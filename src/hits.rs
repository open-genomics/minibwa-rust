//! hit 生成与处理，逐行对照 C 版 `map-algo.c`。

use crate::dna::{hash64, hash_str};
use crate::l2bit::L2b;
use crate::seed::Anchor;

pub const PARENT_UNSET: i32 = -1;
pub const PARENT_TMP_PRI: i32 = -2;

/// 比对附加信息（对应 `mb_extra_t`，M5 碱基比对后填充）。
#[derive(Debug, Clone, Default)]
pub struct Extra {
    pub dp_score: i32,
    pub dp_max0: i32,
    pub dp_max: i32,
    pub dp_max2: i32,
    pub n_ambi: u32,
    pub cigar: Vec<u32>,
    pub cs: Option<String>,
}

/// 一条 hit（对应 `mb_hit_t`）。
#[derive(Debug, Clone, Default)]
pub struct Hit {
    pub tid: i64,
    pub ts: i64,
    pub te: i64,
    pub id: i32,
    pub cnt: i32,
    pub score: i32,
    pub score0: i32,
    pub a_off: i32, // C 版字段名 as（锚点在数组中的起点）
    pub qs: i32,
    pub qe: i32,
    pub parent: i32,
    pub n_sub: i32,
    pub subsc: i32,
    pub mlen: i32,
    pub blen: i32,
    pub mapq: i32,
    pub hash: u32,
    pub rev: bool,
    pub proper_pair: bool,
    pub sam_pri: bool,
    pub flt: bool,
    pub inv: bool,
    pub split: u8,
    pub split_inv: bool,
    pub rescued: bool,
    pub frac_high: u8,
    pub seed_ratio: u8,
    pub p: Option<Box<Extra>>,
}

/// 计算链的 mlen/blen（对应 `mb_cal_mblen`）。
pub fn cal_mblen(a: &[Anchor]) -> (i32, i32) {
    if a.is_empty() {
        return (0, 0);
    }
    let mut mlen = a[0].len as i64;
    let mut blen = a[0].len as i64;
    for i in 1..a.len() {
        let span = a[i].len as i64;
        let tl = a[i].tpos - a[i - 1].tpos;
        let ql = a[i].qpos as i64 - a[i - 1].qpos as i64;
        blen += tl.max(ql);
        mlen += if tl > span && ql > span {
            span
        } else if tl < ql {
            tl
        } else {
            ql
        };
    }
    (mlen as i32, blen as i32)
}

/// 由链的首尾锚点设置 hit 坐标（对应 `mb_hit_set_coor`）。
fn hit_set_coor(r: &mut Hit, qlen: i32, a: &[Anchor]) {
    let k = r.a_off as usize;
    let ak0 = &a[k];
    let ak1 = &a[k + r.cnt as usize - 1];
    r.tid = (ak0.sid >> 1) as i64;
    r.rev = ak0.sid & 1 != 0;
    r.ts = ak0.tpos + 1 - ak0.len as i64;
    r.te = ak1.tpos + 1;
    if !r.rev {
        r.qs = ak0.qpos + 1 - ak0.len;
        r.qe = ak1.qpos + 1;
    } else {
        r.qs = qlen - (ak1.qpos + 1);
        r.qe = qlen - (ak0.qpos + 1 - ak0.len);
    }
    (r.mlen, r.blen) = cal_mblen(&a[k..k + r.cnt as usize]);
}

/// 高覆盖度区间统计（对应 `mb_cal_high_cov`）。
pub fn cal_high_cov(sai: &[crate::bwt::Sai], max_occ: i32) -> i32 {
    let mut b: Vec<u64> = sai
        .iter()
        .filter(|s| s.size > max_occ as u64)
        .map(|s| s.info)
        .collect();
    if b.is_empty() {
        return 0;
    }
    b.sort_unstable();
    let mut hi_st = (b[0] >> 32) as i64;
    let mut hi_en = (b[0] & 0xffff_ffff) as i64;
    let mut hi_cov = 0i64;
    for &x in &b[1..] {
        let st = (x >> 32) as i64;
        let en = (x & 0xffff_ffff) as i64;
        if st > hi_en {
            hi_cov += hi_en - hi_st;
            hi_st = st;
            hi_en = en;
        } else if hi_en < en {
            hi_en = en;
        }
    }
    hi_cov += hi_en - hi_st;
    hi_cov as i32
}

/// 链 → hit（对应 `mb_gen_hit`）。u[i] = score<<32 | 锚点数。
pub fn gen_hit(hash_seed: u32, qlen: i32, _l2b: &L2b, u: &[u64], a: &[Anchor]) -> Vec<Hit> {
    let n_u = u.len();
    if n_u == 0 {
        return Vec::new();
    }
    // z: (排序键, 锚点偏移<<32 | 锚点数)
    let mut z: Vec<(u64, u64)> = Vec::with_capacity(n_u);
    let mut k = 0usize;
    for &ui in u {
        let h = hash64(hash64(a[k].tpos as u64).wrapping_add(hash64(a[k].qpos as i64 as u64)) ^ hash_seed as u64)
            as u32;
        z.push((ui ^ h as u64, (k as u64) << 32 | (ui & 0xffff_ffff)));
        k += (ui & 0xffff_ffff) as usize;
    }
    // C：基数排序（稳定）后反转 → 大分在前
    z.sort_by_key(|&(x, _)| x);
    z.reverse();

    let mut r: Vec<Hit> = Vec::with_capacity(n_u);
    for (i, &(zx, zy)) in z.iter().enumerate() {
        let mut ri = Hit {
            id: i as i32,
            parent: PARENT_UNSET,
            score: (zx >> 32) as i32,
            score0: (zx >> 32) as i32,
            hash: zx as u32,
            cnt: zy as i32,
            a_off: (zy >> 32) as i32,
            ..Default::default()
        };
        hit_set_coor(&mut ri, qlen, a);
        r.push(ri);
    }
    r
}

/// 设置 parent/次级关系（对应 `mb_set_parent`）。
#[allow(clippy::needless_range_loop)]
pub fn set_parent(mask_level: f32, mask_len: i32, r: &mut [Hit], sub_diff: i32, hard_mask_level: bool) {
    let n = r.len();
    if n == 0 {
        return;
    }
    for i in 0..n {
        r[i].id = i as i32;
    }
    let mut w: Vec<usize> = vec![0];
    r[0].parent = 0;
    for i in 1..n {
        let (si, ei) = (r[i].qs, r[i].qe);
        let mut uncov_len = 0i32;
        let mut n_par = 0;

        if !hard_mask_level {
            let mut cov: Vec<u64> = Vec::new();
            for &j in &w {
                let (mut sj, mut ej) = (r[j].qs, r[j].qe);
                if ej <= si || sj >= ei {
                    continue;
                }
                if sj < si {
                    sj = si;
                }
                if ej > ei {
                    ej = ei;
                }
                cov.push((sj as u64) << 32 | ej as u64);
            }
            if cov.is_empty() {
                // 无覆盖：直接成为 primary
                w.push(i);
                r[i].parent = i as i32;
                r[i].n_sub = 0;
                continue;
            }
            cov.sort_unstable();
            let mut x = si;
            for &c in &cov {
                let st = (c >> 32) as i32;
                let en = (c & 0xffff_ffff) as i32;
                if st > x {
                    uncov_len += st - x;
                }
                if en > x {
                    x = en;
                }
            }
            if ei > x {
                uncov_len += ei - x;
            }
        }

        // 找重叠最大的 parent 候选
        let mut max_ol = 0i32;
        let mut max_j = None;
        for &j in &w {
            let (sj, ej) = (r[j].qs, r[j].qe);
            let ol = if ej <= si || sj >= ei {
                0
            } else {
                ej.min(ei) - sj.max(si)
            };
            if max_ol < ol {
                max_ol = ol;
                max_j = Some(j);
            }
        }
        if let Some(mj) = max_j {
            n_par += update_sub(&mut r[..], i, mj, mask_level, mask_len, sub_diff, uncov_len);
            if n_par == 0 {
                for &j in &w {
                    n_par += update_sub(&mut r[..], i, j, mask_level, mask_len, sub_diff, uncov_len);
                    if n_par > 0 {
                        break;
                    }
                }
            }
        }
        if n_par == 0 {
            w.push(i);
            r[i].parent = i as i32;
            r[i].n_sub = 0;
        }
    }
}

fn update_sub(r: &mut [Hit], i: usize, p: usize, mask_level: f32, mask_len: i32, sub_diff: i32, uncov_len: i32) -> i32 {
    let (si, ei) = (r[i].qs, r[i].qe);
    let (sj, ej) = (r[p].qs, r[p].qe);
    if ej <= si || sj >= ei {
        return 0;
    }
    let min = (ej - sj).min(ei - si);
    let max = (ej - sj).max(ei - si);
    let ol = ej.min(ei) - sj.max(si);
    if (ol as f64) / min as f64 - (uncov_len as f64) / max as f64 > mask_level as f64 && uncov_len <= mask_len {
        let mut cnt_sub = 0;
        let sci = r[i].score;
        let parent_of_p = r[p].parent;
        r[i].parent = parent_of_p;
        if r[p].subsc < sci {
            r[p].subsc = sci;
        }
        if r[p].p.is_some() && r[i].p.is_some()
            && (r[p].tid != r[i].tid || r[p].ts != r[i].ts || r[p].te != r[i].te || ol != min)
        {
            let sci2 = r[i].p.as_ref().unwrap().dp_max;
            let dp_max_p = r[p].p.as_mut().unwrap();
            if dp_max_p.dp_max2 < sci2 {
                dp_max_p.dp_max2 = sci2;
            }
            if dp_max_p.dp_max - r[i].p.as_ref().unwrap().dp_max <= sub_diff {
                cnt_sub = 1;
            }
        }
        if cnt_sub > 0 {
            r[p].n_sub += 1;
        }
        1
    } else {
        0
    }
}

/// 选择 SAM primary（对应 `mb_set_sam_pri`）。返回 primary 下标。
#[allow(clippy::needless_range_loop)]
pub fn set_sam_pri(r: &mut [Hit], is_primary5: bool) -> i32 {
    let n = r.len();
    if n == 0 {
        return -1;
    }
    let mut n_pri = 0;
    let mut min_i = -1i32;
    let mut min_qs = -1i32;
    let mut first_i = -1i32;
    for i in 0..n {
        r[i].sam_pri = false;
        if r[i].id != r[i].parent {
            continue;
        }
        n_pri += 1;
        if n_pri == 1 {
            first_i = i as i32;
        }
        if min_qs < 0 || r[i].qs < min_qs {
            min_i = i as i32;
            min_qs = r[i].qs;
        }
    }
    assert!(n_pri > 0);
    let new_pri = if is_primary5 { min_i } else { first_i };
    r[new_pri as usize].sam_pri = true;
    new_pri
}

/// 保留次优 hit（对应 `mb_select_sub`）。
#[allow(clippy::needless_range_loop)]
pub fn select_sub(pri_ratio: f32, min_diff: i32, best_n: i32, r: &mut Vec<Hit>) {
    if pri_ratio <= 0.0 || r.is_empty() {
        return;
    }
    let n = r.len();
    let mut keep = vec![false; n];
    let mut n_2nd = 0i32;
    for i in 0..n {
        let p = r[i].parent as usize;
        if p == i || r[i].inv {
            keep[i] = true;
        } else if (r[i].score as f32 >= r[p].score as f32 * pri_ratio
            || r[i].score + min_diff >= r[p].score)
            && n_2nd < best_n
            && !(r[i].qs == r[p].qs
                && r[i].qe == r[p].qe
                && r[i].tid == r[p].tid
                && r[i].ts == r[p].ts
                && r[i].te == r[p].te)
            {
                keep[i] = true;
                n_2nd += 1;
            }
    }
    let mut k = 0usize;
    for i in 0..n {
        if keep[i] {
            r.swap(k, i);
            k += 1;
        }
    }
    r.truncate(k);
    if k != n {
        sync_hits(r);
    }
}

/// 重排后同步 id/parent（对应 `mb_sync_hits`）。
pub fn sync_hits(regs: &mut [Hit]) {
    let n = regs.len();
    if n == 0 {
        return;
    }
    let mut max_id = -1i32;
    for r in regs.iter() {
        if r.id > max_id {
            max_id = r.id;
        }
    }
    let mut tmp = vec![-1i32; max_id as usize + 1];
    for (i, r) in regs.iter().enumerate() {
        if r.id >= 0 {
            tmp[r.id as usize] = i as i32;
        }
    }
    for (i, r) in regs.iter_mut().enumerate() {
        r.id = i as i32;
        if r.parent == PARENT_TMP_PRI {
            r.parent = i as i32;
        } else if r.parent >= 0 && tmp[r.parent as usize] >= 0 {
            r.parent = tmp[r.parent as usize];
        } else {
            r.parent = PARENT_UNSET;
        }
    }
    set_sam_pri(regs, false);
}

/// 计算 mapping quality（对应 `mb_set_mapq`）。
#[allow(clippy::needless_range_loop)]
pub fn set_mapq(qlen: i32, regs: &mut [Hit], min_chain_sc: i32, match_sc: i32, is_sr: bool, max_sr_len: i32) {
    const MAPQ_COEF_LEN: f64 = 50.0;
    const MAPQ_COEF_FAC: f64 = 3.0; // bwa-mem 的 bug 值（应为 ln(50)），保持一致
    const Q_COEF: f64 = 40.0;
    if regs.is_empty() {
        return;
    }
    for i in 0..regs.len() {
        let inv = regs[i].inv;
        let is_pri = regs[i].parent == regs[i].id;
        if inv {
            regs[i].mapq = 0;
            continue;
        }
        if !is_pri {
            regs[i].mapq = 0;
            continue;
        }
        let r = &regs[i];
        let score = r.score;
        let score0 = r.score0;
        let n_sub = r.n_sub;
        let seed_ratio = r.seed_ratio;
        let mlen = r.mlen;
        let blen = r.blen;
        let subsc_raw = r.subsc;
        let p_info = r.p.as_ref().map(|p| (p.dp_max, p.dp_max2));

        let pen_chn = if score as f64 > qlen as f64 * 0.1 {
            1.0
        } else {
            10.0 * score as f64 / qlen as f64
        };
        let subsc = if subsc_raw > min_chain_sc { subsc_raw } else { min_chain_sc };
        let mut mapq;
        if let Some((dp_max, dp_max2)) = p_info {
            if dp_max2 > 0 && dp_max > 0 {
                let identity = mlen as f64 / blen as f64;
                // bwa-mem 短读公式
                let mut x = if (blen as f64) < MAPQ_COEF_LEN {
                    1.0
                } else {
                    MAPQ_COEF_FAC / (blen as f64).ln()
                };
                x *= identity * identity;
                let mapq_sr = (6.02 * x * x * (dp_max - dp_max2) as f64 / match_sc as f64 + 0.499) as i32;
                // minimap2 长读公式
                let mut x = dp_max2 as f64 / dp_max as f64;
                if subsc > score0 {
                    x *= subsc as f64 / score0 as f64;
                }
                let mapq_lr = (pen_chn * identity * Q_COEF * (1.0 - x * x) * (dp_max as f64 / match_sc as f64).ln()) as i32;
                mapq = if is_sr {
                    mapq_sr
                } else if max_sr_len < 0 {
                    mapq_lr
                } else if qlen < max_sr_len {
                    mapq_sr
                } else {
                    (mapq_lr as f64
                        - (mapq_lr - mapq_sr) as f64 * 2.0f64.powf(1.0 - qlen as f64 / max_sr_len as f64)
                        + 0.499) as i32
                };
            } else {
                let x = subsc as f64 / score0 as f64;
                let identity = mlen as f64 / blen as f64;
                mapq = (pen_chn * identity * Q_COEF * (1.0 - x) * (dp_max as f64 / match_sc as f64).ln()) as i32;
            }
        } else {
            let x = subsc as f64 / score0 as f64;
            mapq = (pen_chn * Q_COEF * (1.0 - x) * (score as f64).ln()) as i32;
        }
        mapq -= (4.343 * ((n_sub + 1) as f64).ln() + 0.499) as i32;
        if mapq < 0 {
            mapq = 0;
        }
        if seed_ratio < 50 {
            mapq = (mapq as f64 * seed_ratio as f64 * seed_ratio as f64 / 2500.0) as i32;
        }
        mapq = mapq.min(60);
        if let Some((dp_max, dp_max2)) = p_info {
            if dp_max > dp_max2 && mapq == 0 {
                mapq = 1;
            }
        }
        regs[i].mapq = mapq;
    }
    set_inv_mapq(regs);
}

/// 倒位 hit 的 mapQ 修正（对应 `mb_set_inv_mapq`）。
fn set_inv_mapq(regs: &mut [Hit]) {
    let n = regs.len();
    if n < 3 {
        return;
    }
    let has_inv = regs.iter().any(|r| r.inv);
    if !has_inv {
        return;
    }
    let mut aux: Vec<(u64, usize)> = Vec::new();
    for (i, r) in regs.iter().enumerate() {
        if r.parent == r.id || r.parent < 0 {
            aux.push(((r.tid as u64) << 32 | r.ts as u64, i));
        }
    }
    aux.sort_by_key(|&(x, _)| x);
    for i in 1..aux.len() - 1 {
        let idx = aux[i].1;
        if regs[idx].inv {
            let l = regs[aux[i - 1].1].mapq;
            let r = regs[aux[i + 1].1].mapq;
            regs[idx].mapq = l.min(r);
        }
    }
}

/// qname 哈希（对应 mb_map_sai 开头的 hash 计算）。
pub fn query_hash(qname: Option<&str>, qlen: i64, seed: i32) -> u32 {
    let mut hash: u32 = qname.map(|s| hash_str(s.as_bytes())).unwrap_or(0);
    hash ^= hash64(qlen as u64).wrapping_add(hash64(seed as u64)) as u32;
    hash64(hash as u64) as u32
}
