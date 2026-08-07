//! 单序列映射编排，逐行对照 C 版 `mb_map` / `mb_map_sai`。

use crate::bwt::Sai;
use crate::chain::lchain_dp;
use crate::dna::NT4;
use crate::hits::{
    cal_high_cov, gen_hit, query_hash, select_sub, set_mapq, set_parent, Hit,
};
use crate::index::Idx;
use crate::opt::Opt;
use crate::seed::{anchor_sort, anchors, seed_intv};

/// 单序列映射结果。
#[derive(Default)]
pub struct MapOut {
    pub hits: Vec<Hit>,
    pub hi_cov: i32,
    pub seed_ratio: f64,
}

/// 对应 `mb_map`（非 meth）：自适应参数 + 两轮 SMEM + mb_map_sai。
pub fn map_one(opt: &Opt, idx: &Idx, qname: &str, seq0: &[u8]) -> MapOut {
    let qlen = seq0.len() as i32;
    let opt_adap = opt.adap(qlen);
    let seq: Vec<u8> = seq0.iter().map(|&b| NT4[b as usize]).collect();
    let u = seed_intv(&idx.bwt, &seq, opt.min_len, opt.max_sub_occ);
    map_sai(&opt_adap, idx, qlen, qname, &seq, u)
}

/// 对应 `mb_map_sai`。当前实现 NO_ALN 路径（链级），M5 接入碱基比对。
fn map_sai(opt: &Opt, idx: &Idx, qlen: i32, qname: &str, seq: &[u8], mut u: Vec<Sai>) -> MapOut {
    const MIN_RECHAIN_LEN: i64 = 1000;
    const MIN_RECHAIN_RATIO: f64 = 0.1;
    let mut out = MapOut::default();
    if u.is_empty() {
        return out;
    }
    let hash = query_hash(Some(qname), qlen as i64, opt.seed);
    let hi_cov = cal_high_cov(&u, opt.max_occ);
    let is_sr = opt.is_sr_mode(qlen);
    out.hi_cov = hi_cov;
    out.seed_ratio = 1.0; // 非 meth 模式恒为 1.0

    // 锚点
    let chn_pen_gap = opt.chain_gap_scale * 0.01 * opt.min_len as f32;
    let v = anchors(&idx.l2b, &idx.bwt, &mut u, qlen, opt.max_occ);

    // 初次链化
    let (mut a, mut w) = match lchain_dp(
        &idx.l2b,
        opt.max_gap,
        opt.max_gap,
        opt.bw,
        opt.max_chain_skip,
        opt.max_chain_iter,
        opt.min_chain_score,
        chn_pen_gap,
        v,
    ) {
        Some(x) => x,
        None => return out,
    };
    let n_hit = w.len();

    // 长读重链化
    if opt.bw_long > opt.bw * 2 && !is_sr && n_hit > 0 {
        // 严格大于：保留第一个最高分链（与 C 一致）
        let mut best = 0usize;
        for i in 1..n_hit {
            if (w[i] >> 32) > (w[best] >> 32) {
                best = i;
            }
        }
        let mut as_off = 0usize;
        for &wi in w.iter().take(best) {
            as_off += (wi & 0xffff_ffff) as usize;
        }
        let st = a[as_off].qpos + 1 - a[as_off].len;
        let en = a[as_off + (w[best] & 0xffff_ffff) as usize - 1].qpos + 1;
        let span = (en - st) as i64;
        if qlen as i64 - span > MIN_RECHAIN_LEN && span as f64 > qlen as f64 * MIN_RECHAIN_RATIO {
            anchor_sort(&idx.l2b, &mut a);
            if let Some((a2, w2)) = lchain_dp(
                &idx.l2b,
                opt.max_gap,
                opt.max_gap,
                opt.bw_long,
                opt.max_chain_skip,
                opt.max_chain_iter,
                opt.min_chain_score,
                chn_pen_gap,
                std::mem::take(&mut a),
            ) {
                a = a2;
                w = w2;
            }
        }
    }

    // hit 生成与次优筛选
    let mut hit = gen_hit(hash, qlen, &idx.l2b, &w, &a);
    let sub_diff = (opt.a + opt.b).max(opt.q + opt.e);
    set_parent(opt.mask_level, opt.mask_len, &mut hit, sub_diff, false);
    select_sub(opt.pri_ratio, opt.min_len * 2, opt.best_n, &mut hit);

    // 碱基级比对
    if opt.flag & crate::opt::flag::NO_ALN == 0 {
        crate::align::align_skeleton(opt, idx, qlen, seq, &mut hit, &mut a);
        let sub_diff = (opt.a + opt.b).max(opt.q + opt.e);
        crate::hits::set_parent(opt.mask_level, opt.mask_len, &mut hit, sub_diff, false);
        crate::hits::select_sub(opt.pri_ratio, opt.min_len * 2, opt.best_n, &mut hit);
        crate::hits::set_sam_pri(&mut hit, opt.flag & crate::opt::flag::PRIMARY5 != 0);
    }

    let frac_high = (255.0 * hi_cov as f64 / qlen as f64) as u8;
    let mut seed_ratio8 = (255.0 * out.seed_ratio + 0.499) as u8;
    if seed_ratio8 == 0 {
        seed_ratio8 = 1;
    }
    for r in hit.iter_mut() {
        r.frac_high = frac_high;
        r.seed_ratio = seed_ratio8;
    }
    set_mapq(qlen, &mut hit, opt.min_chain_score, opt.a, is_sr, opt.max_sr_len);

    out.hits = hit;
    out
}
