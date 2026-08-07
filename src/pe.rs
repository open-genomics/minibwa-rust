//! 配对端算法，逐行对照 C 版 `pe.c`。
//!
//! 包含三个层次：
//! - `pestat`：从唯一定位的 read pair 估计 insert size 分布（4 个方向）
//! - `pair_hits` / `matesw*`：proper pair 配对与 mate rescue
//! - `pair`：配对决策（选配对/选单端）与 PE mapQ 计算

use crate::dna::NT4;
use crate::hits::{set_mapq, set_parent, set_sam_pri, Extra, Hit};
use crate::ksw2;
use crate::l2bit::L2b;
use crate::opt::{flag as of, Opt};

/// 配对 insert size 分布（对应 `mb_pestat_t`）。
#[derive(Debug, Clone, Copy)]
pub struct Pestat {
    pub lo: i32,
    pub hi: i32,
    pub avg: f64,
    pub std: f64,
    pub failed: bool,
}

impl Default for Pestat {
    fn default() -> Self {
        Pestat {
            lo: 0,
            hi: 0,
            avg: 0.0,
            std: 0.0,
            failed: true,
        }
    }
}

/// 预定义的 PE 分布（对应 C 版 `s->pes[1] = ...`：仅 FR 方向可用）。
pub fn pes_predefined(opt: &Opt) -> [Pestat; 4] {
    let mut pes = [Pestat::default(); 4];
    pes[1].failed = false;
    pes[1].avg = opt.pe_avg as f64;
    pes[1].std = opt.pe_std as f64;
    pes[1].lo = opt.pe_lo;
    pes[1].hi = opt.pe_hi;
    pes
}

/// 配对辅助信息（对应 `mb_pairaux_t`）。
#[derive(Debug, Default)]
struct PairAux {
    score: i32,
    sub_sc: i32,
    n_sub: i32,
    n_pp: i32,
    i: [i32; 2],
}

/// 两端方向的相对方向与距离（对应 `mb_insert_dir`）。
fn insert_dir(h0: &Hit, h1: &Hit) -> (usize, i64) {
    let p0 = if h0.rev { h0.te } else { h0.ts };
    let p1 = if h1.rev { h1.te } else { h1.ts };
    let dist = if p0 > p1 { p0 - p1 } else { p1 - p0 };
    let dir = (((h0.rev as i32) << 1 | h1.rev as i32) ^ if p0 < p1 { 0 } else { 3 }) as usize;
    (dir, dist)
}

/// 配对打分（对应 `mb_pair_score`）。
fn pair_score(h0: &Hit, h1: &Hit, pes: &[Pestat; 4], match_sc: i32) -> f64 {
    const MB_SQRT1_2: f64 = std::f64::consts::FRAC_1_SQRT_2;
    let (dir, dist) = insert_dir(h0, h1);
    if !pes[dir].failed && dist >= pes[dir].lo as i64 && dist <= pes[dir].hi as i64 {
        let ns = (dist as f64 - pes[dir].avg) / pes[dir].std; // normalized score
        let dp0 = h0.p.as_ref().map_or(h0.score as f64, |p| p.dp_max as f64);
        let dp1 = h1.p.as_ref().map_or(h1.score as f64, |p| p.dp_max as f64);
        return dp0 + dp1 + 0.721 * (2.0 * libm::erfc(ns.abs() * MB_SQRT1_2)).ln() * match_sc as f64;
        // .721 = 1/log(4)
    }
    -1.0
}

/// 选取唯一命中（对应 `mb_select_unique_se`）：仅一个 primary 且 mapQ≥10。
fn select_unique_se(hit: &[Hit]) -> Option<&Hit> {
    let mut n_pri = 0;
    let mut mapq = 0;
    let mut k = None;
    for (j, h) in hit.iter().enumerate() {
        if h.id == h.parent {
            n_pri += 1;
            mapq = h.mapq;
            k = Some(j);
        }
    }
    if n_pri == 1 && mapq >= 10 {
        k.map(|j| &hit[j])
    } else {
        None
    }
}

/// 全部 primary hit 的覆盖分数和（对应 `mb_hit_sum_score`）。
fn hit_sum_score(n_hit: usize, hit: &[Hit]) -> i32 {
    let mut a: Vec<u64> = Vec::new();
    for (i, h) in hit.iter().enumerate().take(n_hit) {
        if h.id == h.parent && h.p.is_some() {
            a.push((h.qs as u64) << 32 | i as u64);
        }
    }
    let n_pri = a.len();
    if n_pri == 0 {
        return 0;
    }
    if n_pri == 1 {
        return hit[(a[0] & 0xffff_ffff) as usize].p.as_ref().unwrap().dp_max;
    }
    a.sort_unstable();
    let mut sc = 0i32;
    let mut qe = 0i32;
    for &x in &a {
        let h = &hit[(x & 0xffff_ffff) as usize];
        if h.qe <= qe {
            continue;
        }
        let dp_max = h.p.as_ref().unwrap().dp_max;
        if h.qs < qe {
            sc += ((h.qe - qe) as f64 / (h.qe - h.qs) as f64 * dp_max as f64 + 0.499) as i32;
        } else {
            sc += dp_max;
        }
        qe = h.qe;
    }
    sc
}

/// 估计 insert size 分布（对应 `mb_pestat`）。
pub fn pestat(
    opt: &Opt,
    n_frag: usize,
    seg_off: &[usize],
    seg_cnt: &[usize],
    seg_hits: &[Vec<Hit>],
    pes: &mut [Pestat; 4],
) {
    const MIN_DIR_CNT: usize = 20;
    const MIN_DIR_RATIO: f64 = 0.05;
    const OUTLIER_BOUND: f64 = 2.0;
    const MAPPING_BOUND: f64 = 3.0;
    const MAX_STDDEV: f64 = 4.0;
    let mut is: [Vec<i64>; 4] = Default::default();
    for p in pes.iter_mut() {
        *p = Pestat::default();
    }
    for i in 0..n_frag {
        if seg_cnt[i] != 2 {
            continue;
        }
        let off = seg_off[i];
        let (r0, r1) = match (select_unique_se(&seg_hits[off]), select_unique_se(&seg_hits[off + 1])) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };
        if r0.tid != r1.tid {
            continue; // not on the same contig
        }
        let (dir, dist) = insert_dir(r0, r1);
        if dist < opt.max_pe_ins as i64 {
            is[dir].push(dist);
        }
    }
    let max = is.iter().map(Vec::len).max().unwrap_or(0);
    for d in 0..4 {
        let q = &mut is[d];
        if q.len() < MIN_DIR_CNT || q.len() < (max as f64 * MIN_DIR_RATIO) as usize {
            pes[d].failed = true;
            continue;
        }
        q.sort_unstable();
        pes[d].failed = false;
        let n = q.len();
        let p25 = q[((0.25 * n as f64) + 0.499) as usize];
        let p50 = q[((0.50 * n as f64) + 0.499) as usize];
        let p75 = q[((0.75 * n as f64) + 0.499) as usize];
        let r = &mut pes[d];
        r.lo = ((p25 as f64 - OUTLIER_BOUND * (p75 - p25) as f64 + 0.499) as i32).max(1);
        r.hi = (p75 as f64 + OUTLIER_BOUND * (p75 - p25) as f64 + 0.499) as i32;
        let mut x = 0usize;
        let mut avg = 0.0f64;
        for &v in q.iter() {
            if v >= r.lo as i64 && v <= r.hi as i64 {
                avg += v as f64;
                x += 1;
            }
        }
        r.avg = avg / x as f64;
        let mut std = 0.0f64;
        for &v in q.iter() {
            if v >= r.lo as i64 && v <= r.hi as i64 {
                std += (v as f64 - r.avg) * (v as f64 - r.avg);
            }
        }
        r.std = (std / x as f64).sqrt();
        r.lo = (p25 as f64 - MAPPING_BOUND * (p75 - p25) as f64 + 0.499) as i32;
        r.hi = (p75 as f64 + MAPPING_BOUND * (p75 - p25) as f64 + 0.499) as i32;
        if (r.lo as f64) > r.avg - MAX_STDDEV * r.std {
            r.lo = (r.avg - MAX_STDDEV * r.std + 0.499) as i32;
        }
        if (r.hi as f64) < r.avg + MAX_STDDEV * r.std {
            r.hi = (r.avg + MAX_STDDEV * r.std + 0.499) as i32;
        }
        if r.lo < 1 {
            r.lo = 1;
        }
        let _ = p50;
    }
}

/// 扫描 proper pairs（对应 `mb_pair_hits`）。返回最优配对信息。
#[allow(clippy::needless_range_loop)]
fn pair_hits(opt: &Opt, l2b: &L2b, n_hit: &[i32; 2], hit: &mut [Vec<Hit>; 2], pes: &[Pestat; 4]) -> PairAux {
    let mut ret = PairAux {
        score: -1,
        sub_sc: -1,
        n_sub: 0,
        n_pp: 0,
        i: [-1, -1],
    };
    if n_hit[0] == 0 || n_hit[1] == 0 {
        return ret;
    }
    // pa: (x = 参考坐标, y = i<<2 | rev<<1 | r)
    let mut pa: Vec<(i64, u64)> = Vec::with_capacity((n_hit[0] + n_hit[1]) as usize);
    for r in 0..2usize {
        for i in 0..n_hit[r] as usize {
            hit[r][i].proper_pair = false;
            let h = &hit[r][i];
            let x = l2b.ctg[h.tid as usize].off as i64 + if h.rev { h.te } else { h.ts };
            let y = ((i as u64) << 2) | ((h.rev as u64) << 1) | (r as u64);
            pa.push((x, y));
        }
    }
    pa.sort_unstable_by_key(|p| p.0);
    let mut y = [-1i32; 4];
    let mut pp: Vec<(u64, u64)> = Vec::new(); // pp: proper pairs
    let mut marks: Vec<(usize, usize)> = Vec::new(); // (read 端, hit 下标)
    for i in 0..pa.len() {
        let (pi_x, pi_y) = pa[i];
        let hi_r = (pi_y & 1) as usize;
        let hi_i = (pi_y >> 2) as usize;
        let hi_rev = ((pi_y >> 1) & 1) as usize;
        let hi_tid = hit[hi_r][hi_i].tid;
        for r in 0..2usize {
            let dir = (r << 1) | hi_rev;
            if pes[dir].failed {
                continue; // invalid orientation
            }
            let which = (r << 1) | (hi_r ^ 1);
            if y[which] < 0 {
                continue; // no previous hit
            }
            let mut k = y[which] as isize;
            while k >= 0 {
                let (pk_x, pk_y) = pa[k as usize];
                if (pk_y & 3) != which as u64 {
                    k -= 1;
                    continue;
                }
                let hk_r = (pk_y & 1) as usize;
                let hk_i = (pk_y >> 2) as usize;
                if hi_tid != hit[hk_r][hk_i].tid {
                    break;
                }
                let dist = pi_x - pk_x;
                if dist > pes[dir].hi as i64 {
                    break;
                }
                if dist < pes[dir].lo as i64 {
                    k -= 1;
                    continue;
                }
                marks.push((hi_r, hi_i));
                marks.push((hk_r, hk_i));
                let hk_hash = hit[hk_r][hk_i].hash;
                let hi_hash = hit[hi_r][hi_i].hash;
                let mut s = pair_score(&hit[hk_r][hk_i], &hit[hi_r][hi_i], pes, opt.a);
                if s < 0.0 {
                    s = 0.0;
                }
                let qy = if (pk_y & 1) == 0 {
                    ((pk_y >> 2) << 32) | (pi_y >> 2)
                } else {
                    ((pi_y >> 2) << 32) | (pk_y >> 2)
                }; // 高 32 位：read[0] 的 hit 下标
                let qx = ((s + 0.499) as u64) << 32 | (hk_hash ^ hi_hash) as u64;
                pp.push((qx, qy));
                k -= 1;
            }
        }
        y[(pi_y & 3) as usize] = i as i32;
    }
    for (r, i) in marks {
        hit[r][i].proper_pair = true;
    }
    ret.n_pp = pp.len() as i32;
    if ret.n_pp > 0 {
        let tmp = (opt.a + opt.b).max(opt.q + opt.e);
        let mut max = 0u64;
        let mut max2 = 0u64;
        for &(qx, qy) in &pp {
            if qx > max {
                max2 = max;
                max = qx;
                ret.i[0] = (qy >> 32) as i32;
                ret.i[1] = (qy & 0xffff_ffff) as i32;
            } else if qx > max2 {
                max2 = qx;
            }
        }
        assert!(ret.i[0] < n_hit[0] && ret.i[1] < n_hit[1]);
        ret.score = (max >> 32) as i32;
        ret.sub_sc = (max2 >> 32) as i32;
        for &(qx, _) in &pp {
            if (qx >> 32) >= (ret.score - tmp) as u64 {
                ret.n_sub += 1;
            }
        }
    }
    ret
}

/// 无 gap 的线性查找（对应 `mb_ungap`）：kmer 的 Hough 变换找候选对角线。
/// 返回 (max, max_i, n_good, n_kmer)。
#[allow(clippy::needless_range_loop)]
fn ungap(qlen: usize, qseq: &[u8], tlen: usize, tseq: &[u8], kmer: i32) -> (i32, i32, i32, i32) {
    let cap = 1usize << (kmer * 2);
    let mask = cap - 1;
    let mut h = vec![0u16; cap];
    let mut a = vec![0i32; tlen];
    let mut n_kmer = 0i32;
    let mut max_i = -1i32;
    if qlen >= u16::MAX as usize {
        return (0, -1, 0, 0);
    }
    let mut l = 0usize;
    let mut x = 0u32;
    for i in 0..qlen {
        let c = qseq[i];
        if c < 4 {
            x = (x << 2 | c as u32) & mask as u32;
            l += 1;
            if l >= kmer as usize {
                if h[x as usize] == 0 {
                    n_kmer += 1; // 不同 k-mer 数
                }
                h[x as usize] = i as u16;
            }
        } else {
            x = 0;
            l = 0;
        }
    }
    let mut l = 0usize;
    let mut x = 0u32;
    for i in 0..tlen {
        let c = tseq[i];
        if c < 4 {
            x = (x << 2 | c as u32) & mask as u32;
            l += 1;
            if l >= kmer as usize {
                let hv = h[x as usize] as i32;
                if hv > 0 && (i as i32) >= hv {
                    a[(i as i32 - hv) as usize] += 1;
                }
            }
        } else {
            x = 0;
            l = 0;
        }
    }
    let mut max = 0i32;
    for i in 0..tlen {
        if max < a[i] {
            max = a[i];
            max_i = i as i32;
        }
    }
    let mut n_good = 0i32;
    for i in 0..tlen {
        if a[i] > max >> 1 {
            n_good += 1;
        }
    }
    (max, max_i, n_good, n_kmer)
}

/// mate 的碱基级比对（对应 `mb_matesw_align`）。qseq/tseq 为 nt4 编码，原地反转。
fn matesw_align(
    opt: &Opt,
    qlen: i32,
    qseq: &mut [u8],
    tlen: i32,
    tseq: &mut [u8],
    h: &mut Hit,
    min_sc: i32,
    ez: &mut ksw2::Extz,
) {
    let mut mat = [0i8; 25];
    let max_sc = qlen.min(tlen);
    let b_mm = (opt.b + opt.a - 1) / opt.a;
    let b_ts = (opt.b_ts + opt.a - 1) / opt.a;
    let b_ambi = (opt.b_ambi + opt.a - 1) / opt.a;
    let gapo = (opt.q + opt.a - 1) / opt.a;
    let gape = (opt.e + opt.a - 1) / opt.a;
    *h = Hit::default();
    if max_sc >= 32767 {
        return;
    }
    ksw2::gen_nt4_mat(&mut mat, 1, b_mm as i8, b_ts as i8, b_ambi as i8, 0);
    let rst = ksw2::ll_core(qseq, tseq, 5, &mat, gapo, gape, opt.min_len);
    if rst.score >= opt.min_dp_max && rst.score >= min_sc {
        // min_sc 已除以 opt->a
        let te = rst.te + 1;
        let qe = rst.qe + 1;
        let mut qseq2 = qseq[..qe as usize].to_vec();
        let mut tseq2 = tseq[..te as usize].to_vec();
        qseq2.reverse(); // mb_seq_rev
        tseq2.reverse();
        ksw2::gen_nt4_mat(&mut mat, opt.a as i8, opt.b as i8, opt.b_ts as i8, opt.b_ambi as i8, 0);
        let ksw_flag = ksw2::flag::EXTZ_ONLY | ksw2::flag::RIGHT | ksw2::flag::REV_CIGAR;
        ksw2::ksw_extd2(
            qe,
            &qseq2,
            te,
            &tseq2,
            5,
            &mat,
            opt.q as i8,
            opt.e as i8,
            opt.q as i8,
            opt.e as i8,
            opt.bw,
            opt.zdrop,
            opt.end_bonus,
            ksw_flag,
            ez,
        );
        if !ez.cigar.is_empty() && ez.max >= opt.min_dp_max * opt.a {
            let mut p = Box::new(Extra::default());
            p.cigar = std::mem::take(&mut ez.cigar);
            h.p = Some(p);
            h.rescued = true;
            h.qe = qe;
            h.te = te as i64;
            h.ts = te as i64
                - if ez.reach_end {
                    ez.mqe_t as i64 + 1
                } else {
                    ez.max_t as i64 + 1
                };
            h.qs = qe - if ez.reach_end { qe } else { ez.max_q + 1 };
            {
                let p = h.p.as_mut().unwrap();
                p.dp_max = ez.max;
                p.dp_score = ez.max;
                p.dp_max2 = ((ez.max as f64 / rst.score as f64) * rst.score2 as f64 + 0.499) as i32;
                if p.dp_max2 < 0 {
                    p.dp_max2 = 0;
                }
            }
            h.score = rst.score;
            h.score0 = rst.score;
            h.subsc = rst.score2;
            h.cnt = 0;
            h.a_off = -1;
            h.parent = crate::hits::PARENT_UNSET;
            // C 版在调用 mb_update_extra 前已把 qseq/tseq 反转恢复为原方向
            qseq2.reverse();
            tseq2.reverse();
            crate::align::update_extra(
                h,
                &qseq2[h.qs as usize..],
                &tseq2[h.ts as usize..],
                &mat,
                opt.q as i8,
                opt.e as i8,
                opt.flag,
                false,
            );
        }
    }
}

/// rescue 核心（对应 `mb_matesw_core`）：对候选 hit 在允许方向区间内做 SW。
/// 成功时向 h1 追加 rescued hit，返回其下标。
fn matesw_core(
    opt: &Opt,
    l2b: &L2b,
    pes: &[Pestat; 4],
    h0: &Hit,
    r0: usize,
    len: usize,
    seq: &[Vec<u8>; 2], // seq[0]=正链 nt4, seq[1]=反链 nt4
    h1: &mut Vec<Hit>,
    min_sc: i32,
    ez: &mut ksw2::Extz,
) -> Option<usize> {
    let mut skip = [false; 4];
    for dir in 0..4 {
        skip[dir] = pes[dir].failed;
    }
    if skip[0] && skip[1] && skip[2] && skip[3] {
        return None; // no need to perform SW
    }
    let pos5 = if h0.rev { h0.te } else { h0.ts };
    let ctg_len = l2b.ctg[h0.tid as usize].len as i64;
    for dir in 0..4 {
        if skip[dir] {
            continue;
        }
        let is_rev = ((dir >> 1 != dir & 1) as usize) ^ (h0.rev as usize) != 0; // 是否反转互补 mate
        let is_larger = if dir >> 1 != dir & 1 {
            (dir >> 1) ^ (is_rev as usize) != 0 // FR/RF：mate 坐标更大？
        } else {
            (dir >> 1) ^ r0 ^ (h0.rev as usize) != 0 // FF/RR
        };
        let mut ts = (if is_larger {
            pos5 + pes[dir].lo as i64
        } else {
            pos5 - pes[dir].hi as i64
        }) - if is_rev { len as i64 } else { 0 };
        let mut te = (if is_larger {
            pos5 + pes[dir].hi as i64
        } else {
            pos5 - pes[dir].lo as i64
        }) + if !is_rev { len as i64 } else { 0 };
        if ts < 0 {
            ts = 0;
        }
        if te > ctg_len {
            te = ctg_len;
        }
        if te - ts > len as i64 {
            let mut ts2 = ts;
            let mut te2 = te;
            let mut ref_seq = vec![0u8; (te - ts) as usize];
            l2b.getseq(h0.tid, ts, te, &mut ref_seq);
            let (max_ug, max_i, n_good, n_kmer) = ungap(len, &seq[is_rev as usize], (te - ts) as usize, &ref_seq, 7);
            if max_ug >= 10 && max_ug >= (len as i32) >> 1 && n_good == 1 {
                ts2 = ts + max_i as i64 - (len as i64) / 2;
                te2 = ts2 + (len as i64) * 2;
                if ts2 < ts {
                    ts2 = ts;
                }
                if te2 > te {
                    te2 = te;
                }
            }
            if max_ug >= 10 || (max_ug as f64) >= n_kmer as f64 * 0.33 {
                let mut qseq = seq[is_rev as usize][..len].to_vec();
                let mut tseq = ref_seq[(ts2 - ts) as usize..(te2 - ts) as usize].to_vec();
                let mut ht = Hit::default();
                matesw_align(
                    opt,
                    len as i32,
                    &mut qseq,
                    (te2 - ts2) as i32,
                    &mut tseq,
                    &mut ht,
                    min_sc,
                    ez,
                );
                if ht.p.is_some() {
                    ht.tid = h0.tid;
                    ht.ts += ts2;
                    ht.te += ts2;
                    ht.rev = is_rev;
                    if is_rev {
                        let qt = ht.qs;
                        ht.qs = len as i32 - ht.qe;
                        ht.qe = len as i32 - qt;
                    }
                    h1.push(ht);
                    return Some(h1.len() - 1);
                }
            }
        }
    }
    None
}

/// mate rescue（对应 `mb_matesw`）。返回新增的 hit 数。
fn matesw(
    opt: &Opt,
    l2b: &L2b,
    n_hit: &mut [i32; 2],
    hit: &mut [Vec<Hit>; 2],
    pes: &[Pestat; 4],
    paux0: &PairAux,
    qlen: &[i32; 2],
    qseq: &[Vec<u8>; 2], // ASCII
) -> i32 {
    if opt.max_rescue == 0 {
        return 0;
    }
    // 预计算 rescue 候选数
    let mut n_res = 0i32;
    for r in 0..2 {
        let h0 = &hit[r];
        let mut m = 0i32;
        for h in h0.iter().take(n_hit[r] as usize) {
            if m >= opt.max_rescue {
                break;
            }
            if !h.proper_pair
                && h.p.as_ref().unwrap().dp_max >= h0[0].p.as_ref().unwrap().dp_max - opt.pen_unpair * opt.a
            {
                m += 1;
            }
        }
        n_res += m;
    }
    if n_res == 0 {
        return 0;
    }
    // 收集 rescue 候选（必须与上面的循环一致）
    let mut a: Vec<(u64, u64)> = Vec::with_capacity(n_res as usize);
    for r in 0..2 {
        let h0 = &hit[r];
        let mut m = 0i32;
        for (i, h) in h0.iter().enumerate().take(n_hit[r] as usize) {
            if m >= opt.max_rescue {
                break;
            }
            if !h.proper_pair
                && h.p.as_ref().unwrap().dp_max >= h0[0].p.as_ref().unwrap().dp_max - opt.pen_unpair * opt.a
            {
                a.push((
                    (h.p.as_ref().unwrap().dp_max as u64) << 32 | h.hash as u64,
                    (i as u64) << 1 | r as u64,
                ));
                m += 1;
            }
        }
    }
    a.sort_unstable_by_key(|p| p.0);
    // 准备比对序列（nt4 正链 / 反链）
    let mut qs: [[Vec<u8>; 2]; 2] = Default::default();
    for r in 0..2 {
        for i in 0..qlen[r] as usize {
            let c = NT4[qseq[r][i] as usize];
            qs[r][0].push(c);
            qs[r][1].push(if c < 4 { 3 - c } else { 4 });
        }
        qs[r][1].reverse();
    }
    let mut min_sc = [0i32; 2];
    for r in 0..2 {
        min_sc[r] = hit_sum_score(n_hit[r] as usize, &hit[r]) / opt.a - opt.pen_unpair;
    }
    // rescue（从高分到低分）
    let mut max = [paux0.score; 2];
    let mut max2 = [paux0.sub_sc; 2];
    let mut skip = [false; 2];
    let mut ez = ksw2::Extz::default();
    for i in (0..a.len()).rev() {
        let r = (a[i].1 & 1) as usize;
        let j = (a[i].1 >> 1) as usize;
        if skip[r] {
            continue;
        }
        let h0 = hit[r][j].clone();
        let found = matesw_core(
            opt,
            l2b,
            pes,
            &h0,
            r,
            qlen[1 - r] as usize,
            &qs[1 - r],
            &mut hit[1 - r],
            min_sc[1 - r],
            &mut ez,
        );
        if found.is_some() {
            let h1 = hit[1 - r].last().unwrap();
            let sc = pair_score(&h0, h1, pes, opt.a) as i32;
            if sc > max[r] {
                max2[r] = max[r];
                max[r] = sc;
            } else if sc > max2[r] {
                max2[r] = sc;
            }
            if max[r] == max2[r] {
                skip[r] = true;
            }
        }
    }
    let n_add = (hit[0].len() - n_hit[0] as usize + hit[1].len() - n_hit[1] as usize) as i32;
    for r in 0..2 {
        n_hit[r] = hit[r].len() as i32;
    }
    n_add
}

/// 按分数降序排序 hit（对应 `mb_hit_sort` 的排序语义）。
fn hit_sort(r: &mut [Hit]) {
    if r.len() <= 1 {
        return;
    }
    r.sort_by(|a, b| {
        let sa = a.p.as_ref().map_or(a.score, |p| p.dp_max);
        let sb = b.p.as_ref().map_or(b.score, |p| p.dp_max);
        sb.cmp(&sa).then(a.hash.cmp(&b.hash))
    });
}

/// 高覆盖分数同步（对应 `mb_sync_high_cov`）。
fn sync_high_cov(h: &mut [Hit]) {
    let mut max_frac = 0u8;
    for x in h.iter() {
        if x.frac_high > max_frac {
            max_frac = x.frac_high;
        }
    }
    for x in h.iter_mut() {
        x.frac_high = max_frac;
    }
}

/// 配对（对应 `mb_pair`）。修改 hit 的 mapq/parent/proper_pair/sam_pri 等。
pub fn pair(
    opt: &Opt,
    l2b: &L2b,
    n_hit: &mut [i32; 2],
    hit: &mut [Vec<Hit>; 2],
    pes: &[Pestat; 4],
    qlen: &[i32; 2],
    qseq: &[Vec<u8>; 2], // ASCII
) {
    const PE_BONUS: i32 = 4;
    if n_hit[0] == 0 && n_hit[1] == 0 {
        return;
    }
    let seed_ratio = [
        if n_hit[0] > 0 { hit[0][0].seed_ratio } else { 255 },
        if n_hit[1] > 0 { hit[1][0].seed_ratio } else { 255 },
    ];
    let min_seed_ratio = seed_ratio[0].min(seed_ratio[1]);
    let mut paux = pair_hits(opt, l2b, n_hit, hit, pes);
    let do_matesw = !(paux.n_pp > 0 && paux.score == paux.sub_sc); // 两个并列最优则跳过 rescue
    if do_matesw && opt.max_rescue > 0 {
        let sub_diff = (opt.a + opt.b).max(opt.q + opt.e);
        if matesw(opt, l2b, n_hit, hit, pes, &paux, qlen, qseq) > 0 {
            #[allow(clippy::needless_range_loop)]
            for r in 0..2 {
                for h in hit[r].iter_mut() {
                    if !h.rescued {
                        h.n_sub = 0;
                        h.subsc = 0;
                        if let Some(p) = h.p.as_mut() {
                            p.dp_max2 = 0;
                        }
                    }
                }
                hit_sort(&mut hit[r]);
                set_parent(opt.mask_level, opt.mask_len, &mut hit[r], sub_diff, false);
                set_mapq(
                    qlen[r],
                    &mut hit[r],
                    opt.min_chain_score,
                    opt.a,
                    opt.is_sr_mode(qlen[r]),
                    opt.max_sr_len,
                );
            }
            paux = pair_hits(opt, l2b, n_hit, hit, pes); // 有 rescued hit 时重新配对
        }
    }
    if paux.n_pp == 0 {
        // end_pairing：无 proper pair，恢复 sam_pri
        set_sam_pri(&mut hit[0], opt.flag & of::PRIMARY5 != 0);
        set_sam_pri(&mut hit[1], opt.flag & of::PRIMARY5 != 0);
        return;
    }
    let mut dp_max_se = [0i32; 2];
    let mut dp_max_se2 = [0i32; 2];
    for r in 0..2 {
        for h in hit[r].iter() {
            let dm = h.p.as_ref().map_or(0, |p| p.dp_max);
            if dp_max_se[r] < dm {
                dp_max_se2[r] = dp_max_se[r];
                dp_max_se[r] = dm;
            } else if dp_max_se2[r] < dm {
                dp_max_se2[r] = dm;
            }
        }
    }
    let score_se = dp_max_se[0] + dp_max_se[1];
    let score_se2 = dp_max_se2[0] + dp_max_se2[1];
    sync_high_cov(&mut hit[0]);
    sync_high_cov(&mut hit[1]);
    for h in hit.iter_mut() {
        // 清 0x2：下面会重新评估
        for h in h.iter_mut() {
            h.proper_pair = false;
        }
    }
    if paux.score >= score_se - opt.pen_unpair * opt.a {
        // 选配对
        let hi = [paux.i[0] as usize, paux.i[1] as usize];
        let identity =
            (hit[0][hi[0]].mlen + hit[1][hi[1]].mlen) as f64 / (hit[0][hi[0]].blen + hit[1][hi[1]].blen) as f64;
        let mut score2 = paux.sub_sc;
        if (hit[0][hi[0]].id != hit[0][hi[0]].parent || hit[1][hi[1]].id != hit[1][hi[1]].parent)
            && score2 < score_se - opt.pen_unpair * opt.a
        {
            score2 = score_se - opt.pen_unpair * opt.a;
        }
        let mut diff = paux.score - score2;
        if diff > paux.score + PE_BONUS * opt.a - score_se2 {
            // inspired by dragmap
            diff = paux.score + PE_BONUS * opt.a - score_se2;
        }
        let mut mapq_pe = (6.02 * identity * identity * diff as f64 / opt.a as f64
            - 4.343 * ((paux.n_sub + 1) as f64).ln()
            + 0.499) as i32;
        if mapq_pe < 0 {
            mapq_pe = 0;
        }
        mapq_pe = (mapq_pe as f64
            * (1.0 - 0.5 * (hit[0][hi[0]].frac_high as f64 / 255.0 + hit[1][hi[1]].frac_high as f64 / 255.0))
            + 0.499) as i32;
        if min_seed_ratio < 50 {
            mapq_pe = (mapq_pe as f64 * min_seed_ratio as f64 * min_seed_ratio as f64 / 2500.0) as i32;
        }
        if mapq_pe > 60 {
            mapq_pe = 60;
        }
        if mapq_pe <= 0 && paux.score > score2 {
            mapq_pe = 1;
        }
        for r in 0..2 {
            let q = hi[r];
            if n_hit[r] == 1 {
                // 单 hit：取 max(mapq_se, mapq_pe)
                hit[r][q].mapq = hit[r][q].mapq.max(mapq_pe);
            } else {
                // 多 hit：封顶并惩罚小 mapq_se
                hit[r][q].mapq = if hit[r][q].mapq > mapq_pe {
                    mapq_pe
                } else {
                    (0.2 * hit[r][q].mapq as f64 + 0.8 * mapq_pe as f64 + 0.499) as i32
                };
            }
            if hit[r][q].id != hit[r][q].parent {
                // 提升配对 hit 为 primary 并更新 parent
                let p_id = hit[r][q].parent;
                for i in 0..hit[r].len() {
                    if hit[r][i].parent == p_id {
                        hit[r][i].parent = hit[r][q].id;
                    }
                }
                hit[r][p_id as usize].mapq = 0;
            }
            for i in 0..hit[r].len() {
                // 处理其他嵌合 hit
                if i == q {
                    continue;
                }
                if hit[r][i].id == hit[r][i].parent {
                    let (si, ei) = (hit[r][i].qs, hit[r][i].qe);
                    let (sq, eq) = (hit[r][q].qs, hit[r][q].qe);
                    let ol = if ei <= sq || si >= eq {
                        0
                    } else {
                        ei.min(eq) - si.max(sq)
                    };
                    if (ol as f64) > opt.mask_level as f64 * (ei - si) as f64 {
                        for j in 0..hit[r].len() {
                            if hit[r][j].parent == hit[r][i].id {
                                hit[r][j].parent = hit[r][q].id;
                            }
                        }
                        hit[r][i].mapq = 0;
                    }
                }
            }
            for h in hit[r].iter_mut() {
                h.sam_pri = false;
                h.proper_pair = false;
                if h.id == h.parent {
                    h.proper_pair = true;
                }
            }
            hit[r][q].sam_pri = true;
        }
        if opt.flag & of::PRIMARY5 != 0 {
            let pri = [set_sam_pri(&mut hit[0], true), set_sam_pri(&mut hit[1], true)];
            if pri[0] != paux.i[0] || pri[1] != paux.i[1] {
                // sam_pri 变化则清除 0x2
                for h in hit.iter_mut() {
                    for h in h.iter_mut() {
                        h.proper_pair = false;
                    }
                }
            }
        }
    } else {
        // 选单端
        let diff = score_se - opt.pen_unpair * opt.a - paux.score;
        let mapq_pe = 6 * diff / opt.a;
        for h in hit.iter_mut() {
            for h in h.iter_mut() {
                h.mapq = h.mapq.min(mapq_pe);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pestat_default_is_failed() {
        let p = Pestat::default();
        assert!(p.failed);
    }

    #[test]
    fn ll_core_basic() {
        // query/target 完全一致（nt4）
        let q: Vec<u8> = b"ACGTACGTACGT".iter().map(|&b| NT4[b as usize]).collect();
        let t: Vec<u8> = b"TTACGTACGTACGT".iter().map(|&b| NT4[b as usize]).collect();
        let mut mat = [0i8; 25];
        ksw2::gen_nt4_mat(&mut mat, 1, 4, 4, 1, 0);
        let r = crate::ksw2::ll_core(&q, &t, 5, &mat, 4, 1, 0);
        assert_eq!(r.score, 12); // 12 个匹配
        assert_eq!(r.te, 13); // target 2..14，末尾坐标 13
        assert_eq!(r.qe, 11);
    }
}
