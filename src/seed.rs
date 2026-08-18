//! 种子：两轮 SMEM 与 anchor 生成，逐行对照 C 版 `seed.c`。

use crate::bwt::{Bwt, Sai};
use crate::l2bit::L2b;

/// 锚点（对应 `mb_anchor_t`）。
pub const SEED_LONG_JOIN: u32 = 0x1;
pub const SEED_IGNORE: u32 = 0x2;

#[derive(Debug, Clone, Copy, Default)]
pub struct Anchor {
    /// tid<<1 | rev
    pub sid: i32,
    /// 锚点长度
    pub len: i32,
    /// 查询坐标（锚点最后一个碱基；rev 时已翻转）
    pub qpos: i32,
    /// 目标坐标（contig 内）
    pub tpos: i64,
    pub flt: bool,
    pub flag: u32,
}

/// 两轮 SMEM（对应 `mb_seed_intv` 的单序列版）。
///
/// pass 1：标准 SMEM（min_occ=1）；
/// pass 2：对长度 ≥ 2*min_len 且 occ ≤ max_sub_occ 的 SMEM，
/// 在其区间内以更小的 min_len、更高的 min_occ 找子 SMEM。
pub fn seed_intv(bwt: &Bwt, seq: &[u8], min_len: i32, max_sub_occ: i32) -> Vec<Sai> {
    let len = seq.len() as i64;
    let mut v: Vec<Sai> = Vec::new();
    let mut p = Sai::default();

    let mut x = 0i64;
    while {
        x = bwt.smem(seq, x, min_len as i64, 1, &mut p);
        if p.size > 0 {
            v.push(p);
        }
        x < len
    } {}

    let n_a0 = v.len();
    for i in 0..n_a0 {
        let st = (v[i].info >> 32) as i64;
        let en = (v[i].info & 0xffff_ffff) as i64;
        let parent_size = v[i].size;
        if en - st < (min_len as i64) * 2 || parent_size > max_sub_occ as u64 {
            continue;
        }
        let sub_min_len = ((en - st) / 2).max(min_len as i64);
        let mut x = st;
        while {
            x = bwt.smem(seq, x, sub_min_len, parent_size as i64 + 1, &mut p);
            if p.size > parent_size {
                v.push(p);
            }
            x < en
        } {}
    }
    v
}

/// 区间排序去重（对应 `mb_seed_sort_dedup`）：
/// 按 x[0] 升序、size 降序、info 升序排序，去除三字段全同的重复。
pub fn seed_sort_dedup(u: &mut Vec<Sai>) {
    if u.len() <= 1 {
        return;
    }
    u.sort_by(|a, b| a.x[0].cmp(&b.x[0]).then(b.size.cmp(&a.size)).then(a.info.cmp(&b.info)));
    u.dedup_by(|a, b| a.x[0] == b.x[0] && a.size == b.size && a.info == b.info);
}

/// 种子区间 → 锚点（对应 `mb_anchor` 的非 meth 路径）。
///
/// 对每个区间解析全部（或均匀采样至 max_occ 个）SA 值，
/// 换算为 contig 坐标；按拼接坐标排序后去除被包含的短锚点。
pub fn anchors(l2b: &L2b, bwt: &Bwt, u: &mut Vec<Sai>, qlen: i32, max_occ: i32) -> Vec<Anchor> {
    let mut v: Vec<Anchor> = Vec::new();
    if u.is_empty() {
        return v;
    }
    seed_sort_dedup(u);

    // 按 (x[0], size) 分组处理（对应 aux[]）
    let mut i0 = 0usize;
    let n = u.len();
    let mut i = 1usize;
    while i0 < n {
        while i < n && u[i].x[0] == u[i0].x[0] && u[i].size == u[i0].size {
            i += 1;
        }
        let q = u[i0];
        // 收集 BWT 行号：size <= max_occ 全取，否则均匀采样
        let mut rows: Vec<u64> = Vec::new();
        if q.size <= max_occ as u64 {
            for j in 0..q.size {
                rows.push(q.x[0] + j);
            }
        } else {
            let mut j = 0u64;
            let mut taken = 0i32;
            while j < q.size && taken < max_occ {
                let mut step = (q.size - j) / (max_occ - taken) as u64;
                if step < 1 {
                    step = 1;
                }
                rows.push(q.x[0] + j);
                j += step;
                taken += 1;
            }
        }
        // 每个 SA 值 × 组内每个区间 → 锚点
        for &row in &rows {
            let sa = bwt.sa(row);
            for u_j in &u[i0..i] {
                let qs = (u_j.info >> 32) as i32;
                let qe = (u_j.info & 0xffff_ffff) as i32;
                let len = qe - qs;
                let mut rev = false;
                let (tid, cst) = l2b.intv2cid(sa, sa + len as u64, &mut rev);
                if tid < 0 {
                    continue;
                }
                let ctg = &l2b.ctg[tid as usize];
                let rev_i = i32::from(rev);
                v.push(Anchor {
                    sid: ((tid as i32) << 1) | rev_i,
                    len,
                    qpos: if rev { qlen - 1 - qs } else { qs + len - 1 },
                    tpos: (ctg.off * 2 + ctg.len * rev as u64) as i64 + cst + (len - 1) as i64,
                    flt: false,
                    flag: 0,
                });
            }
        }
        i0 = i;
        i += 1;
    }

    // 按拼接坐标 tpos 稳定排序（C 版 LSD 基数排序是稳定的）
    v.sort_by_key(|a| a.tpos);
    // 调整为 contig 内坐标
    for a in v.iter_mut() {
        let ctg = &l2b.ctg[(a.sid >> 1) as usize];
        a.tpos -= (ctg.off * 2 + ctg.len * (a.sid & 1) as u64) as i64;
    }
    anchor_dedup(&mut v);
    v
}

/// 去除被更长锚点完全包含的短锚点（对应 `mb_anchor_dedup`，要求按 tpos 有序）。
fn anchor_dedup(v: &mut Vec<Anchor>) {
    const MAX_BACK: i64 = 100;
    let n = v.len() as i64;
    for i in 1..n {
        let ai = v[i as usize];
        let tsi = ai.tpos + 1 - ai.len as i64;
        let qsi = ai.qpos + 1 - ai.len;
        let mut k = 0i64;
        let mut j = i - 1;
        while j >= 0 && k < MAX_BACK {
            let aj = v[j as usize];
            if aj.sid != ai.sid {
                break;
            }
            if aj.tpos < tsi {
                break;
            }
            let tsj = aj.tpos + 1 - aj.len as i64;
            let qsj = aj.qpos + 1 - aj.len;
            if tsj >= tsi {
                // j 被 i 包含
                if tsj - tsi == (qsj - qsi) as i64 {
                    v[j as usize].flt = true;
                }
            } else if ai.tpos == aj.tpos {
                // i 被 j 包含
                if ai.qpos == aj.qpos {
                    v[i as usize].flt = true;
                }
            }
            j -= 1;
            k += 1;
        }
    }
    v.retain(|a| !a.flt);
}

/// 锚点重排序：临时加上拼接偏移，按 tpos 稳定排序，再减回（对应 `mb_anchor_sort`）。
pub fn anchor_sort(l2b: &L2b, a: &mut [Anchor]) {
    if a.len() <= 1 {
        return;
    }
    for x in a.iter_mut() {
        let ctg = &l2b.ctg[(x.sid >> 1) as usize];
        x.tpos += (ctg.off * 2 + ctg.len * (x.sid & 1) as u64) as i64;
    }
    a.sort_by_key(|x| x.tpos);
    for x in a.iter_mut() {
        let ctg = &l2b.ctg[(x.sid >> 1) as usize];
        x.tpos -= (ctg.off * 2 + ctg.len * (x.sid & 1) as u64) as i64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index;

    /// 构建小型索引（走真实流水线）用于种子测试。
    fn build_toy() -> (L2b, Bwt) {
        let fa = ">c1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
                  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        let l2b = L2b::import(fa.as_bytes(), 11).unwrap();
        let bwt = index::bwt_libsais(&l2b, 4, true, false, 1).unwrap();
        (l2b, bwt)
    }

    #[test]
    fn seed_intv_finds_exact_seeds() {
        let (_, bwt) = build_toy();
        // 查询是参考的前 30bp：应找到覆盖它的 SMEM
        let q: Vec<u8> = (0..30).map(bwt_test_base).collect();
        let mut v = seed_intv(&bwt, &q, 10, 10);
        assert!(!v.is_empty());
        // 至少一个区间覆盖查询起点且足够长
        let best = v.iter().max_by_key(|p| p.size).unwrap();
        let st = best.info >> 32;
        let en = best.info & 0xffff_ffff;
        assert_eq!(st, 0);
        assert!(en >= 20);
        seed_sort_dedup(&mut v);
        // 去重后无三字段全同
        for w in v.windows(2) {
            assert!(w[0].x[0] != w[1].x[0] || w[0].size != w[1].size || w[0].info != w[1].info);
        }
    }

    fn bwt_test_base(i: u64) -> u8 {
        // ACGT 重复
        [0, 1, 2, 3][(i % 4) as usize]
    }

    #[test]
    fn anchors_map_to_expected_positions() {
        let (l2b, bwt) = build_toy();
        let q: Vec<u8> = (0..30).map(bwt_test_base).collect();
        let mut v = seed_intv(&bwt, &q, 10, 10);
        let a = anchors(&l2b, &bwt, &mut v, 30, 250);
        assert!(!a.is_empty());
        // 参考为 ACGT 重复：正向锚点应在 c1 上，tpos 与 qpos 对齐
        let fwd: Vec<&Anchor> = a.iter().filter(|x| x.sid & 1 == 0).collect();
        assert!(!fwd.is_empty());
        for x in &fwd {
            assert_eq!(x.sid >> 1, 0);
            // ACGT 周期为 4：各合法映射位置与查询同相
            assert_eq!((x.tpos - x.qpos as i64).rem_euclid(4), 0);
        }
        // 首个映射应与查询同起点
        assert_eq!(fwd[0].tpos, fwd[0].qpos as i64);
    }
}
