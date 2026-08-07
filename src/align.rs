//! 碱基级比对流水线，逐行对照 C 版 `align.c`。

use crate::dna::log2_approx;
use crate::hits::{cal_mblen, Hit, Extra, PARENT_TMP_PRI};
use crate::index::Idx;
use crate::ksw2::{self, flag as kf, Extz, gen_nt4_mat, CIGAR_MATCH, CIGAR_INS, CIGAR_DEL};
use crate::l2bit::L2b;
use crate::opt::{flag as of, Opt};
use crate::seed::{Anchor, SEED_IGNORE, SEED_LONG_JOIN};

const CIGAR_STR: &[u8; 9] = b"MIDNSHP=X";

fn min_int(a: i32, b: i32) -> i32 { if a < b { a } else { b } }

fn max_bw_from_mm(opt: &Opt, mm: i32) -> i32 {
    let x = mm * (opt.a + opt.b);
    let mut max1 = 0;
    let mut max2 = 0;
    if x >= opt.q + opt.e { max1 = (x - opt.q + opt.e - 1) / opt.e; }
    if x >= opt.q2 + opt.e2 { max2 = (x - opt.q2 + opt.e2 - 1) / opt.e2; }
    max1.max(max2)
}

/// ksw 比对封装（对应 `mb_align_pair`）。
fn align_pair(
    opt: &Opt, qlen: i32, qseq: &[u8], tlen: i32, tseq: &[u8],
    mat: &[i8; 25], w: i32, end_bonus: i32, zdrop: i32, ksw_flag: u32, ez: &mut Extz,
) {
    const MAX_BW_ADJ_LEN: i32 = 100;
    let mut n_mm = -1i32;
    let mut ksw_flag = ksw_flag;
    if (opt.b_ts != 0 && opt.b != opt.b_ts) || opt.flag & of::METH != 0 {
        ksw_flag |= kf::GENERIC_SC;
    }
    // ungapped 快捷路径
    if ksw_flag & kf::EXTZ_ONLY != 0 && tlen >= qlen {
        ksw2::reset_extz(ez);
        ez.score = 0;
        ez.max = 0;
        let mut mx = 0i32;
        for j in 0..qlen {
            ez.score += mat[tseq[j as usize] as usize * 5 + qseq[j as usize] as usize] as i32;
            n_mm += i32::from(tseq[j as usize] > 3 || qseq[j as usize] > 3 || mat[tseq[j as usize] as usize * 5 + qseq[j as usize] as usize] < 0);
            if ez.max < ez.score { ez.max = ez.score; ez.max_q = j; ez.max_t = j; }
            mx = mx.max(ez.score);
        }
        let _ = mx;
        if n_mm <= 2 {
            ez.mqe = ez.score;
            ez.mqe_t = qlen - 1;
            if ez.mqe + end_bonus >= ez.max {
                ez.reach_end = true;
                ksw2::push_cigar_vec(&mut ez.cigar, CIGAR_MATCH, qlen as u32);
                return;
            }
        }
    } else if qlen == tlen && ksw_flag & kf::EXTZ_ONLY == 0 {
        let max_gapped = (qlen - 2) * opt.a - 2 * (opt.q + opt.e);
        ksw2::reset_extz(ez);
        ez.score = 0;
        for j in 0..qlen {
            ez.score += mat[tseq[j as usize] as usize * 5 + qseq[j as usize] as usize] as i32;
            n_mm += i32::from(tseq[j as usize] > 3 || qseq[j as usize] > 3 || mat[tseq[j as usize] as usize * 5 + qseq[j as usize] as usize] < 0);
        }
        if n_mm <= 3 || ez.score > max_gapped {
            ksw2::push_cigar_vec(&mut ez.cigar, CIGAR_MATCH, qlen as u32);
            return;
        }
    }
    // 带宽调整
    let mut w = w;
    if n_mm >= 0 && min_int(qlen, tlen) < MAX_BW_ADJ_LEN {
        let max_bw = max_bw_from_mm(opt, n_mm);
        if w > max_bw + 4 { w = max_bw + 4; }
    }
    // 内存检查 + DP
    if opt.max_sw_mat > 0 && (tlen as i64) * (qlen as i64) > opt.max_sw_mat {
        ksw2::reset_extz(ez);
        ez.zdropped = true;
    } else {
        ksw2::ksw_extd2(qlen, qseq, tlen, tseq, 5, mat, opt.q as i8, opt.e as i8, opt.q2 as i8, opt.e2 as i8,
                         w, zdrop * opt.a, end_bonus, ksw_flag, ez);
    }
}

fn collect_long_gaps(as1: i32, cnt1: i32, a: &[Anchor], min_gap: i64) -> Vec<i32> {
    let mut k = Vec::new();
    for i in 1..cnt1 {
        let gap = (a[(as1 + i) as usize].qpos as i64 - a[(as1 + i - 1) as usize].qpos as i64)
            - (a[(as1 + i) as usize].tpos - a[(as1 + i - 1) as usize].tpos);
        if gap < -min_gap || gap > min_gap { k.push(i); }
    }
    if k.len() <= 1 { return Vec::new(); }
    k
}

fn filter_bad_seeds(as1: i32, cnt1: i32, a: &mut [Anchor], min_gap: i64, diff_thres: i32, max_ext_len: i32, max_ext_cnt: i32) {
    let k = collect_long_gaps(as1, cnt1, a, min_gap);
    if k.is_empty() { return; }
    let n = k.len();
    let mut max = 0i32; let mut max_st = -1i32; let mut max_en = -1i32;
    let mut idx = 0usize;
    loop {
        if idx >= n || (max_en > 0 && idx as i32 >= max_en) {
            if max_en > 0 {
                for i in k[max_st as usize]..k[max_en as usize] {
                    a[(as1 + i) as usize].flag |= SEED_IGNORE;
                }
            }
            max = 0; max_st = -1; max_en = -1;
            if idx >= n { break; }
        }
        let i = k[idx];
        let mut n_ins = 0i32; let mut n_del = 0i32;
        let gap = (a[(as1 + i) as usize].qpos as i64 - a[(as1 + i - 1) as usize].qpos as i64)
            - (a[(as1 + i) as usize].tpos - a[(as1 + i - 1) as usize].tpos);
        if gap > 0 { n_ins += gap as i32; } else { n_del += (-gap) as i32; }
        let qs = a[(as1 + i - 1) as usize].qpos;
        let ts = a[(as1 + i - 1) as usize].tpos;
        let mut max_diff = 0i32; let mut max_diff_l = -1i32;
        for l in (idx + 1)..n.min(idx + 1 + max_ext_cnt as usize) {
            let j = k[l];
            if a[(as1 + j) as usize].qpos - a[(as1 + j) as usize].len - qs > max_ext_len
                || (a[(as1 + j) as usize].tpos - a[(as1 + j) as usize].len as i64 - ts) as i32 > max_ext_len { break; }
            let g = (a[(as1 + j) as usize].qpos as i64 - a[(as1 + j - 1) as usize].qpos as i64)
                - (a[(as1 + j) as usize].tpos - a[(as1 + j - 1) as usize].tpos);
            if g > 0 { n_ins += g as i32; } else { n_del += (-g) as i32; }
            let diff = n_ins + n_del - (n_ins - n_del).abs();
            if max_diff < diff { max_diff = diff; max_diff_l = l as i32; }
        }
        if max_diff > diff_thres && max_diff > max {
            max = max_diff; max_st = idx as i32; max_en = max_diff_l;
        }
        idx += 1;
    }
}

fn filter_bad_seeds_alt(as1: i32, cnt1: i32, a: &mut [Anchor], min_gap: i64, max_ext: i32) {
    let k = collect_long_gaps(as1, cnt1, a, min_gap);
    if k.is_empty() { return; }
    let n = k.len();
    let mut ki = 0usize;
    while ki < n {
        let i = k[ki];
        let mut gap1 = (a[(as1 + i) as usize].qpos as i64 - a[(as1 + i - 1) as usize].qpos as i64)
            - (a[(as1 + i) as usize].tpos - a[(as1 + i - 1) as usize].tpos);
        gap1 = gap1.abs();
        let mut te1 = a[(as1 + i) as usize].tpos;
        let mut qe1 = a[(as1 + i) as usize].qpos;
        let mut left_len = a[(as1 + i) as usize].len;
        let mut l = ki + 1;
        while l < n {
            let j = k[l];
            if a[(as1 + j) as usize].qpos - qe1 > max_ext || (a[(as1 + j) as usize].tpos - te1) as i32 > max_ext { break; }
            let mut gap2 = (a[(as1 + j) as usize].qpos as i64 - a[(as1 + j - 1) as usize].qpos as i64)
                - (a[(as1 + j) as usize].tpos - a[(as1 + j - 1) as usize].tpos);
            let m_t = (a[(as1 + j - 1) as usize].tpos - te1 + left_len as i64) as i32;
            let m_q = a[(as1 + j - 1) as usize].qpos - qe1 + left_len;
            let m = min_int(m_t, m_q);
            gap2 = gap2.abs();
            if m > gap1 as i32 + gap2 as i32 { break; }
            te1 = a[(as1 + j) as usize].tpos;
            qe1 = a[(as1 + j) as usize].qpos;
            left_len = a[(as1 + j) as usize].len;
            gap1 = gap2;
            l += 1;
        }
        if l > ki + 1 {
            let end = k[l - 1];
            for j in k[ki]..end { a[(as1 + j) as usize].flag |= SEED_IGNORE; }
            a[(as1 + end) as usize].flag |= SEED_LONG_JOIN;
        }
        ki = l;
    }
}

fn fix_bad_ends(r: &Hit, a: &[Anchor], bw: i32, min_match: i32) -> (i32, i32) {
    let mut as_ = r.a_off; let mut cnt = r.cnt;
    if r.cnt < 3 { return (as_, cnt); }
    // 左端
    let mut m = a[r.a_off as usize].len;
    let mut l = a[r.a_off as usize].len;
    let mut i = r.a_off + 1;
    while i < r.a_off + r.cnt - 1 {
        let q_span = a[i as usize].len;
        if a[i as usize].flag & SEED_LONG_JOIN != 0 { break; }
        let lr = (a[i as usize].tpos - a[(i - 1) as usize].tpos) as i32;
        let lq = a[i as usize].qpos - a[(i - 1) as usize].qpos;
        let mn = min_int(lr, lq);
        let mx = lr.max(lq);
        if mx - mn > l >> 1 { as_ = i; }
        l += mn;
        m += min_int(mn, q_span);
        if l >= bw << 1 || (m >= min_match && m >= bw) || m >= r.mlen >> 1 { break; }
        i += 1;
    }
    cnt = r.a_off + r.cnt - as_;
    // 右端
    m = a[(r.a_off + r.cnt - 1) as usize].len;
    l = a[(r.a_off + r.cnt - 1) as usize].len;
    i = r.a_off + r.cnt - 2;
    while i > as_ {
        let q_span = a[i as usize].len;
        if a[(i + 1) as usize].flag & SEED_LONG_JOIN != 0 { break; }
        let lr = (a[(i + 1) as usize].tpos - a[i as usize].tpos - a[(i + 1) as usize].len as i64 + a[i as usize].len as i64) as i32;
        let lq = a[(i + 1) as usize].qpos - a[i as usize].qpos - a[(i + 1) as usize].len + a[i as usize].len;
        let mn = min_int(lr, lq);
        let mx = lr.max(lq);
        if mx - mn > l >> 1 { cnt = i + 1 - as_; }
        l += mn;
        m += min_int(mn, q_span);
        if l >= bw << 1 || (m >= min_match && m >= bw) || m >= r.mlen >> 1 { break; }
        i -= 1;
    }
    (as_, cnt)
}

fn max_stretch(r: &Hit, a: &[Anchor]) -> (i32, i32) {
    let mut as_ = r.a_off; let mut cnt = r.cnt;
    if r.cnt < 2 { return (as_, cnt); }
    let mut max_score = -1i32; let mut max_i = -1i32; let mut max_len = 0i32;
    let mut score = a[r.a_off as usize].len; let mut len = 1i32;
    for i in (r.a_off + 1)..(r.a_off + r.cnt) {
        let lr = (a[i as usize].tpos - a[(i - 1) as usize].tpos) as i32;
        let lq = a[i as usize].qpos - a[(i - 1) as usize].qpos;
        if lq == lr {
            score += min_int(lq, a[i as usize].len);
            len += 1;
        } else {
            if score > max_score { max_score = score; max_len = len; max_i = i - len; }
            score = a[i as usize].len; len = 1;
        }
    }
    if score > max_score { max_score = score; max_len = len; max_i = r.a_off + r.cnt - len; }
    as_ = max_i; cnt = max_len;
    (as_, cnt)
}

/// 修正 CIGAR：左对齐 indel、合并零长操作（对应 `mb_fix_cigar`）。
fn fix_cigar(r: &mut Hit, qseq: &[u8], tseq: &[u8]) -> (i32, i32) {
    let p = r.p.as_mut().unwrap();
    let mut qshift = 0i32; let mut tshift = 0i32;
    let mut toff = 0i32; let mut qoff = 0i32;
    let mut to_shrink = false;
    let n = p.cigar.len();
    if n <= 1 { return (0, 0); }
    // indel 左对齐
    for k in 0..n {
        let op = p.cigar[k] & 0xf;
        let len = (p.cigar[k] >> 4) as i32;
        if len == 0 { to_shrink = true; }
        if op == CIGAR_MATCH { toff += len; qoff += len; }
        else if op == CIGAR_INS || op == CIGAR_DEL {
            if k > 0 && k < n - 1 && p.cigar[k - 1] & 0xf == CIGAR_MATCH && p.cigar[k + 1] & 0xf == CIGAR_MATCH {
                let prev_len = (p.cigar[k - 1] >> 4) as i32;
                let mut l = 0;
                if op == CIGAR_INS {
                    while l < prev_len && qseq[(qoff - 1 - l) as usize] != qseq[(qoff + len - 1 - l) as usize] { l += 1; }
                } else {
                    while l < prev_len && tseq[(toff - 1 - l) as usize] != tseq[(toff + len - 1 - l) as usize] { l += 1; }
                }
                if l > 0 {
                    p.cigar[k - 1] -= (l as u32) << 4;
                    p.cigar[k + 1] += (l as u32) << 4;
                    qoff -= l; toff -= l;
                }
                if l == prev_len { to_shrink = true; }
            }
            if op == CIGAR_INS { qoff += len; } else { toff += len; }
        } else if op == 3 { toff += len; } // N_SKIP
    }
    // 修复 5I6D7I 模式
    let mut k = 0;
    while k < n.saturating_sub(2) {
        let op = p.cigar[k] & 0xf;
        if op > 0 && op + (p.cigar[k + 1] & 0xf) == 3 {
            let mut s = [0u32; 3];
            let mut l = k;
            while l < n {
                let o = p.cigar[l] & 0xf;
                if o == CIGAR_INS || o == CIGAR_DEL || p.cigar[l] >> 4 == 0 {
                    s[o as usize] += p.cigar[l] >> 4;
                } else { break; }
                l += 1;
            }
            if s[1] > 0 && s[2] > 0 && l - k > 2 {
                p.cigar[k] = s[1] << 4 | CIGAR_INS;
                p.cigar[k + 1] = s[2] << 4 | CIGAR_DEL;
                for m in (k + 2)..l { p.cigar[m] &= 0xf; }
                to_shrink = true;
            }
            k = l;
        } else { k += 1; }
    }
    if to_shrink {
        let mut l = 0usize;
        for k in 0..n {
            if p.cigar[k] >> 4 != 0 { p.cigar[l] = p.cigar[k]; l += 1; }
        }
        p.cigar.truncate(l);
        let nn = p.cigar.len();
        let mut l = 0usize;
        for k in 0..nn {
            if k == nn - 1 || p.cigar[k] & 0xf != p.cigar[k + 1] & 0xf {
                p.cigar[l] = p.cigar[k]; l += 1;
            } else {
                p.cigar[k + 1] += p.cigar[k] >> 4 << 4;
            }
        }
        p.cigar.truncate(l);
    }
    // 去除前导 I/D
    if !p.cigar.is_empty() {
        let op0 = p.cigar[0] & 0xf;
        if op0 == CIGAR_INS || op0 == CIGAR_DEL {
            let l = (p.cigar[0] >> 4) as i32;
            if op0 == CIGAR_INS {
                if r.rev { r.qe -= l; } else { r.qs += l; }
                qshift = l;
            } else { r.ts += l as i64; tshift = l; }
            p.cigar.remove(0);
        }
    }
    (qshift, tshift)
}

/// 计算 mlen/blen/n_ambi/dp_max，处理 EQX 和 tag（对应 `mb_update_extra`）。
fn update_extra(r: &mut Hit, qseq: &[u8], tseq: &[u8], mat: &[i8; 25], q: i8, e: i8, opt_flag: u64, log_gap: bool) {
    let p = r.p.as_mut().unwrap();
    let (qshift, tshift) = fix_cigar(r, qseq, tseq);
    let qseq = &qseq[qshift as usize..];
    let tseq = &tseq[tshift as usize..];
    r.blen = 0; r.mlen = 0;
    let mut s = 0.0f64; let mut mx = 0.0f64;
    let mut toff = 0i32; let mut qoff = 0i32;
    let p_ref = r.p.as_mut().unwrap();
    let n = p_ref.cigar.len();
    for k in 0..n {
        let op = p_ref.cigar[k] & 0xf;
        let len = (p_ref.cigar[k] >> 4) as i32;
        if op == CIGAR_MATCH {
            let mut n_ambi = 0; let mut n_diff = 0;
            for l in 0..len {
                if (qoff + l) as usize >= qseq.len() || (toff + l) as usize >= tseq.len() { break; }
                let cq = qseq[(qoff + l) as usize]; let ct = tseq[(toff + l) as usize];
                if ct > 3 || cq > 3 { n_ambi += 1; }
                else if ct != cq { n_diff += i32::from(mat[ct as usize * 5 + cq as usize] < 0); }
                s += mat[ct as usize * 5 + cq as usize] as f64;
                if s < 0.0 { s = 0.0; } else if s > mx { mx = s; }
            }
            r.blen += len - n_ambi; r.mlen += len - (n_ambi + n_diff);
            p_ref.n_ambi += n_ambi as u32;
            toff += len; qoff += len;
        } else if op == CIGAR_INS {
            let mut n_ambi = 0;
            for l in 0..len { if (qoff + l) as usize >= qseq.len() { break; } if qseq[(qoff + l) as usize] > 3 { n_ambi += 1; } }
            r.blen += len - n_ambi; p_ref.n_ambi += n_ambi as u32;
            if log_gap { s -= q as f64 + e as f64 * log2_approx(1.0 + len as f32) as f64; }
            else { s -= q as f64 + e as f64; }
            if s < 0.0 { s = 0.0; }
            qoff += len;
        } else if op == CIGAR_DEL {
            let mut n_ambi = 0;
            for l in 0..len { if (toff + l) as usize >= tseq.len() { break; } if tseq[(toff + l) as usize] > 3 { n_ambi += 1; } }
            r.blen += len - n_ambi; p_ref.n_ambi += n_ambi as u32;
            if log_gap { s -= q as f64 + e as f64 * log2_approx(1.0 + len as f32) as f64; }
            else { s -= q as f64 + e as f64; }
            if s < 0.0 { s = 0.0; }
            toff += len;
        } else if op == 3 { toff += len; } // N_SKIP
    }
    p_ref.dp_max0 = (mx + 0.499) as i32;
    p_ref.dp_max = p_ref.dp_max0;
    if opt_flag & of::EQX != 0 { update_cigar_eqx(r, qseq, tseq); }
}

fn update_cigar_eqx(r: &mut Hit, qseq: &[u8], tseq: &[u8]) {
    let p = r.p.as_mut().unwrap();
    let mut new_cigar: Vec<u32> = Vec::new();
    let mut toff = 0i32; let mut qoff = 0i32;
    for k in 0..p.cigar.len() {
        let op = p.cigar[k] & 0xf;
        let len = (p.cigar[k] >> 4) as i32;
        if op == CIGAR_MATCH {
            let mut remaining = len;
            while remaining > 0 {
                let mut l = 0;
                while l < remaining && qseq[(qoff + l) as usize] == tseq[(toff + l) as usize] && qseq[(qoff + l) as usize] < 4 { l += 1; }
                if l > 0 {
                    ksw2::push_cigar_vec(&mut new_cigar, 7, l as u32); // =
                    remaining -= l; toff += l; qoff += l;
                }
                l = 0;
                while l < remaining && !(qseq[(qoff + l) as usize] == tseq[(toff + l) as usize] && qseq[(qoff + l) as usize] < 4) { l += 1; }
                if l > 0 {
                    ksw2::push_cigar_vec(&mut new_cigar, 8, l as u32); // X
                    remaining -= l; toff += l; qoff += l;
                }
            }
        } else if op == CIGAR_INS { qoff += len; ksw2::push_cigar_vec(&mut new_cigar, op, len as u32); }
        else if op == CIGAR_DEL { toff += len; ksw2::push_cigar_vec(&mut new_cigar, op, len as u32); }
        else if op == 3 { toff += len; ksw2::push_cigar_vec(&mut new_cigar, op, len as u32); }
        else { ksw2::push_cigar_vec(&mut new_cigar, op, len as u32); }
    }
    p.cigar = new_cigar;
}

/// 测试 Z-drop 和倒位（对应 `mm_test_zdrop`）。返回 0=无, 1=zdrop, 2=倒位。
fn test_zdrop(opt: &Opt, qseq: &[u8], tseq: &[u8], cigar: &[u32], mat: &[i8; 25], is_sr: bool) -> i32 {
    let mut score = 0i32; let mut mx = i32::MIN; let mut max_i = -1i32; let mut max_j = -1i32;
    let mut i = 0i32; let mut j = 0i32; let mut max_zdrop = 0i32;
    let mut pos = [[-1i32; 2]; 2];
    for &c in cigar {
        let op = c & 0xf; let len = (c >> 4) as i32;
        if op == CIGAR_MATCH {
            for l in 0..len {
                score += mat[tseq[(i + l) as usize] as usize * 5 + qseq[(j + l) as usize] as usize] as i32;
                // update_max_zdrop
                if score < mx {
                    let li = i + l - max_i; let lj = j + l - max_j;
                    let diff = if li > lj { li - lj } else { lj - li };
                    let z = mx - score - diff * opt.e;
                    if z > max_zdrop {
                        max_zdrop = z;
                        pos[0] = [max_i, i + l]; pos[1] = [max_j, j + l];
                    }
                } else { mx = score; max_i = i + l; max_j = j + l; }
            }
            i += len; j += len;
        } else if op == CIGAR_INS || op == CIGAR_DEL || op == 3 {
            score -= opt.q + opt.e * len;
            if op == CIGAR_INS { j += len; } else { i += len; }
            if score < mx {
                let li = i - max_i; let lj = j - max_j;
                let diff = if li > lj { li - lj } else { lj - li };
                let z = mx - score - diff * opt.e;
                if z > max_zdrop { max_zdrop = z; pos[0] = [max_i, i]; pos[1] = [max_j, j]; }
            } else { mx = score; max_i = i; max_j = j; }
        }
    }
    // 倒位检测
    let q_len = pos[1][1] - pos[1][0]; let t_len = pos[0][1] - pos[0][0];
    if !is_sr && max_zdrop > opt.zdrop_inv && q_len < opt.max_gap && t_len < opt.max_gap {
        let qseq2: Vec<u8> = (0..q_len).map(|k| {
            let c = qseq[(pos[1][1] - k - 1) as usize];
            if c >= 4 { 4 } else { 3 - c }
        }).collect();
        let score = ksw2::ll_local(&qseq2, &tseq[pos[0][0] as usize..pos[0][1] as usize], 5, mat, opt.q, opt.e);
        if score >= opt.min_chain_score * opt.a && score >= opt.min_dp_max * opt.a { return 2; }
    }
    if max_zdrop > opt.zdrop { 1 } else { 0 }
}

/// 分割 hit（对应 `mb_split_hit`）。
fn split_hit(r: &mut Hit, r2: &mut Hit, n: i32, qlen: i32, a: &[Anchor], l2b: &L2b) {
    if n <= 0 || n >= r.cnt { return; }
    *r2 = r.clone();
    r2.id = -1; r2.sam_pri = false; r2.p = None;
    r2.split_inv = false; r2.cnt = r.cnt - n;
    r2.score = (r.score as f32 * (r2.cnt as f32 / r.cnt as f32) + 0.499) as i32;
    r2.a_off = r.a_off + n;
    if r.parent == r.id { r2.parent = PARENT_TMP_PRI; }
    hit_set_coor(r2, qlen, a);
    r.cnt -= r2.cnt; r.score -= r2.score;
    hit_set_coor(r, qlen, a);
    r.split |= 1; r2.split |= 2;
}

fn hit_set_coor(r: &mut Hit, qlen: i32, a: &[Anchor]) {
    let k = r.a_off as usize;
    let ak0 = &a[k]; let ak1 = &a[k + r.cnt as usize - 1];
    r.tid = (ak0.sid >> 1) as i64; r.rev = ak0.sid & 1 != 0;
    r.ts = ak0.tpos + 1 - ak0.len as i64; r.te = ak1.tpos + 1;
    if !r.rev { r.qs = ak0.qpos + 1 - ak0.len; r.qe = ak1.qpos + 1; }
    else { r.qs = qlen - (ak1.qpos + 1); r.qe = qlen - (ak0.qpos + 1 - ak0.len); }
    let (mlen, blen) = cal_mblen(&a[k..k + r.cnt as usize]);
    r.mlen = mlen; r.blen = blen;
}

/// 计算 CIGAR 消耗的 query/target 长度。
fn cigar_lens(cigar: &[u32]) -> (i32, i64) {
    let mut q = 0i32; let mut t = 0i64;
    for &c in cigar {
        match c & 0xf {
            0 | 7 | 8 => { let l = (c >> 4) as i32; q += l; t += l as i64; }
            1 => { q += (c >> 4) as i32; }
            2 | 3 => { t += (c >> 4) as i64; }
            _ => {}
        }
    }
    (q, t)
}

/// 核心：单 hit 骨对比对（对应 `mb_align1`）。
fn align1(opt: &Opt, idx: &Idx, qlen: i32, qseq0: &[Vec<u8>], r: &mut Hit, r2: &mut Hit, a: &mut [Anchor]) {
    let is_sr = opt.is_sr_mode(qlen);
    let max_back = if is_sr { 0 } else { 10 };
    r2.cnt = 0;
    if r.cnt == 0 { return; }
    let rev = a[r.a_off as usize].sid & 1 != 0;
    let tid = (a[r.a_off as usize].sid >> 1) as i64;
    let mut mat = [0i8; 25];
    gen_nt4_mat(&mut mat, opt.a as i8, opt.b as i8, opt.b_ts as i8, opt.b_ambi as i8, 0);
    let bw = (opt.bw as f32 * 1.5 + 1.0) as i32;
    let bw_long = if !is_sr {
        let bl = (opt.bw_long as f32 * 1.5 + 1.0) as i32;
        bl.max(bw)
    } else { bw };

    // 确定种子范围
    let (as1, cnt1) = if is_sr {
        max_stretch(r, a)
    } else {
        let (as1, cnt1) = fix_bad_ends(r, a, opt.bw, opt.min_chain_score * 2);
        filter_bad_seeds(as1, cnt1, a, 10, 40, opt.max_gap >> 1, 10);
        filter_bad_seeds_alt(as1, cnt1, a, 30, opt.max_gap >> 1);
        (as1, cnt1)
    };

    // 对齐区域
    let ts = a[as1 as usize].tpos + 1 - a[as1 as usize].len as i64 + min_int(a[as1 as usize].len >> 1, max_back) as i64;
    let qs = a[as1 as usize].qpos + 1 - a[as1 as usize].len + min_int(a[as1 as usize].len >> 1, max_back);
    let te = a[(as1 + cnt1 - 1) as usize].tpos + 1 - min_int(a[(as1 + cnt1 - 1) as usize].len >> 1, max_back) as i64;
    let qe = a[(as1 + cnt1 - 1) as usize].qpos + 1 - min_int(a[(as1 + cnt1 - 1) as usize].len >> 1, max_back);

    let (mut ts0, mut qs0, mut te0, mut qe0);
    if is_sr {
        qs0 = 0; qe0 = qlen;
        let mut l = qs as i64;
        l += if l * (opt.a as i64 + opt.end_bonus as i64) > (opt.q as i64) { (l * opt.a as i64 + opt.end_bonus as i64 - opt.q as i64) / opt.e as i64 } else { 0 };
        ts0 = (ts - l).max(0);
        let l = (qlen - qe) as i64;
        let mut l = l + if l * (opt.a as i64 + opt.end_bonus as i64) > (opt.q as i64) { (l * opt.a as i64 + opt.end_bonus as i64 - opt.q as i64) / opt.e as i64 } else { 0 };
        l = l.min(opt.max_gap as i64);
        te0 = (te + l).min(idx.l2b.ctg[tid as usize].len as i64);
    } else {
        ts0 = (a[r.a_off as usize].tpos + 1 - a[r.a_off as usize].len as i64).max(0);
        qs0 = a[r.a_off as usize].qpos + 1 - a[r.a_off as usize].len;
        let (mut ts1, mut qs1) = (0i64, 0i32);
        if qs > 0 && ts > 0 {
            let l = (qs as i64).min(opt.max_gap as i64);
            qs1 = qs1.max(qs - l as i32);
            let mut l = l + if l * (opt.a as i64) > (opt.q as i64) { (l * opt.a as i64 - opt.q as i64) / opt.e as i64 } else { 0 };
            l = l.min(opt.max_gap as i64).min(ts);
            ts1 = ts1.max(ts - l);
            ts0 = ts0.min(ts1).min(ts);
            qs0 = qs0.min(qs1);
        } else { ts0 = ts; qs0 = qs; let _ = (&mut ts1, &mut qs1); }
        te0 = a[(r.a_off + r.cnt - 1) as usize].tpos + 1;
        qe0 = a[(r.a_off + r.cnt - 1) as usize].qpos + 1;
        let ctg_len = idx.l2b.ctg[tid as usize].len as i64;
        let (mut te1, mut qe1) = (ctg_len, qlen);
        if qe < qlen && te < ctg_len {
            let l = (qlen - qe).min(opt.max_gap);
            qe1 = qe1.min(qe + l);
            qe0 = qe0.max(qe1);
            let mut l = l as i64 + if l as i64 * (opt.a as i64) > (opt.q as i64) { (l as i64 * opt.a as i64 - opt.q as i64) / opt.e as i64 } else { 0 };
            l = l.min(opt.max_gap as i64).min(ctg_len - te);
            te1 = te1.min(te + l);
            te0 = te0.max(te1);
        } else { te0 = te; qe0 = qe; }
    }

    let mut tseq = vec![0u8; (te0 - ts0) as usize];
    let mut ez = Extz::default();
    let mut dropped = false;
    let (mut ts1, mut qs1) = (ts, qs);
    let (mut te1, mut qe1);

    // 初始化 r->p
    if r.p.is_none() {
        r.p = Some(Box::new(Extra::default()));
    }

    // 左延伸
    if qs > 0 && ts > 0 {
        let mut qseg = qseq0[rev as usize][qs0 as usize..qs as usize].to_vec();
        idx.l2b.getseq(tid, ts0, ts, &mut tseq);
        let mut tseg = tseq[..(ts - ts0) as usize].to_vec();
        // 反转后比对（mb_seq_rev 仅反转，不互补）
        qseg.reverse();
        tseg.reverse();
        align_pair(opt, qs - qs0, &qseg, (ts - ts0) as i32, &tseg, &mat, bw, opt.end_bonus,
                   if r.split_inv { opt.zdrop_inv } else { opt.zdrop }, kf::EXTZ_ONLY | kf::RIGHT | kf::REV_CIGAR, &mut ez);
        if !ez.cigar.is_empty() {
            for &c in &ez.cigar { ksw2::push_cigar_vec(&mut r.p.as_mut().unwrap().cigar, c & 0xf, c >> 4); }
            r.p.as_mut().unwrap().dp_score += if ez.reach_end { ez.mqe } else { ez.max };
        }
        ts1 = ts - (if ez.reach_end { ez.mqe_t + 1 } else { ez.max_t + 1 }) as i64;
        qs1 = qs - (if ez.reach_end { qs - qs0 } else { ez.max_q + 1 });
    }
    te1 = ts1; qe1 = qs1;

    // 首个锚点的精确匹配（C 版用 te-ts，不是 te-ts1）
    {
        te1 = a[as1 as usize].tpos + 1 - min_int(a[as1 as usize].len >> 1, max_back) as i64;
        qe1 = a[as1 as usize].qpos + 1 - min_int(a[as1 as usize].len >> 1, max_back);
        let match_len = (te1 - ts) as u32;
        ksw2::push_cigar_vec(&mut r.p.as_mut().unwrap().cigar, CIGAR_MATCH, match_len);
        r.p.as_mut().unwrap().dp_score += opt.a * match_len as i32;
    }
    let (mut ts_cur, mut qs_cur) = (te1, qe1);

    // gap 填充
    for i in 1..cnt1 {
        let ai = &a[(as1 + i) as usize];
        if ai.flag & SEED_IGNORE != 0 && i != cnt1 - 1 { continue; }
        let te1_i = ai.tpos + 1 - min_int(ai.len >> 1, max_back) as i64;
        let qe1_i = ai.qpos + 1 - min_int(ai.len >> 1, max_back);
        if i == cnt1 - 1 || ai.flag & SEED_LONG_JOIN != 0 || (qe1_i - qs_cur >= opt.min_ksw_len && te1_i - ts_cur >= opt.min_ksw_len as i64) {
            let mut bw1 = bw_long;
            let mut d1 = 0i64;
            if ai.len > opt.min_len * 2 {
                d1 = te1_i - (ai.tpos + 1 - ai.len as i64);
                d1 = d1.min((qe1_i - qs_cur) as i64).min(te1_i - ts_cur);
                d1 -= opt.min_len as i64;
                if d1 < opt.min_len as i64 { d1 = 0; }
            }
            let te_use = te1_i - d1; let qe_use = qe1_i - d1 as i32;
            if ai.flag & SEED_LONG_JOIN != 0 {
                bw1 = ((qe_use - qs_cur) as i64).max(te_use - ts_cur) as i32;
            }
            let qseg = qseq0[rev as usize][qs_cur as usize..qe_use as usize].to_vec();
            idx.l2b.getseq(tid, ts_cur, te_use, &mut tseq);
            align_pair(opt, qe_use - qs_cur, &qseg, (te_use - ts_cur) as i32, &tseq[..(te_use - ts_cur) as usize],
                       &mat, bw1, -1, opt.zdrop, kf::APPROX_MAX, &mut ez);
            let zdrop_code = test_zdrop(opt, &qseg, &tseq[..(te_use - ts_cur) as usize], &ez.cigar, &mat, is_sr);
            if zdrop_code != 0 {
                align_pair(opt, qe_use - qs_cur, &qseg, (te_use - ts_cur) as i32, &tseq[..(te_use - ts_cur) as usize],
                           &mat, bw1, -1, if zdrop_code == 2 { opt.zdrop_inv } else { opt.zdrop }, 0, &mut ez);
            }
            if !ez.cigar.is_empty() {
                for &c in &ez.cigar { ksw2::push_cigar_vec(&mut r.p.as_mut().unwrap().cigar, c & 0xf, c >> 4); }
            }
            if ez.zdropped {
                let mut j = (i - 1) - 1;
                while j >= 0 {
                    if a[(as1 + j) as usize].tpos <= ts_cur + ez.max_t as i64 { break; }
                    j -= 1;
                }
                if j < 0 { j = 0; }
                dropped = true;
                r.p.as_mut().unwrap().dp_score += ez.max;
                te1 = ts_cur + (ez.max_t + 1) as i64;
                qe1 = qs_cur + ez.max_q + 1;
                let (mlen, _blen) = cal_mblen(&a[(as1 + j + 1) as usize..(as1 + cnt1) as usize]);
                if mlen >= opt.min_chain_score {
                    split_hit(r, r2, as1 + j + 1 - r.a_off, qlen, a, &idx.l2b);
                    if zdrop_code == 2 { r2.split_inv = true; }
                }
                break;
            } else {
                r.p.as_mut().unwrap().dp_score += ez.score;
            }
            if d1 > 0 {
                ksw2::push_cigar_vec(&mut r.p.as_mut().unwrap().cigar, CIGAR_MATCH, d1 as u32);
                r.p.as_mut().unwrap().dp_score += opt.a * d1 as i32;
                ts_cur = te1_i; qs_cur = qe1_i;
            } else {
                ts_cur = te_use; qs_cur = qe_use;
            }
        }
    }

    // 右延伸
    if !dropped && qs_cur < qe0 && ts_cur < te0 {
        let qseg = qseq0[rev as usize][qs_cur as usize..qe0 as usize].to_vec();
        idx.l2b.getseq(tid, ts_cur, te0, &mut tseq);
        align_pair(opt, qe0 - qs_cur, &qseg, (te0 - ts_cur) as i32, &tseq[..(te0 - ts_cur) as usize],
                   &mat, bw, opt.end_bonus, opt.zdrop, kf::EXTZ_ONLY, &mut ez);
        // 检查 CIGAR 是否覆盖了完整查询长度
        let cigar_q: i32 = ez.cigar.iter().map(|c| {
            let op = c & 0xf; let len = (c >> 4) as i32;
            if op == CIGAR_MATCH || op == CIGAR_INS || op == 7 || op == 8 { len } else { 0 }
        }).sum();
        let force_end = cigar_q == qe0 - qs_cur;
        if !ez.cigar.is_empty() {
            for &c in &ez.cigar { ksw2::push_cigar_vec(&mut r.p.as_mut().unwrap().cigar, c & 0xf, c >> 4); }
            r.p.as_mut().unwrap().dp_score += if ez.reach_end || force_end { ez.mqe.max(ez.score) } else { ez.max };
        }
        te1 = ts_cur + (if ez.reach_end || force_end { (te0 - ts_cur) } else { (ez.max_t + 1) as i64 }) as i64;
        qe1 = qs_cur + (if ez.reach_end || force_end { qe0 - qs_cur } else { ez.max_q + 1 });
    }

    r.ts = ts1; r.te = te1;
    if !rev { r.qs = qs1; r.qe = qe1; }
    else { r.qs = qlen - qe1; r.qe = qlen - qs1; }

    if r.p.is_some() {
        // 同步坐标与 CIGAR：确保 qe-qs == CIGAR query 消耗, te-ts == CIGAR target 消耗
        let (cq, ct) = cigar_lens(&r.p.as_ref().unwrap().cigar);
        let ctg_len = idx.l2b.ctg[tid as usize].len as i64;
        if !rev {
            r.qs = qs1.max(0);
            r.qe = (r.qs + cq).min(qlen);
        } else {
            r.qe = (qlen - qs1.max(0)).min(qlen);
            r.qs = (r.qe - cq).max(0);
        }
        r.ts = ts1.max(0);
        r.te = (r.ts + ct).min(ctg_len);
        let qe_use = if !rev { r.qe } else { qlen - r.qs };
        let qs_use = if !rev { r.qs } else { qlen - r.qe };
        if qs_use >= 0 && qe_use > qs_use && qe_use <= qlen && r.te > r.ts {
            idx.l2b.getseq(tid, r.ts, r.te, &mut tseq);
            let qseq = &qseq0[rev as usize][qs_use as usize..qe_use as usize];
            let tlen_aln = (r.te - r.ts) as usize;
            if tlen_aln <= tseq.len() && qseq.len() == cq as usize {
                update_extra(r, qseq, &tseq[..tlen_aln], &mat, opt.q as i8, opt.e as i8, opt.flag, !is_sr);
            }
        }
    }
}

/// 骨对比对编排（对应 `mb_align_skeleton`）。
pub fn align_skeleton(opt: &Opt, idx: &Idx, qlen: i32, seq: &[u8], hits: &mut Vec<Hit>, a: &mut Vec<Anchor>) {
    let qseq0_0: Vec<u8> = seq.to_vec();
    let qseq0_1: Vec<u8> = {
        let mut v = vec![0u8; qlen as usize];
        for i in 0..qlen as usize {
            v[qlen as usize - 1 - i] = if qseq0_0[i] < 4 { 3 - qseq0_0[i] } else { 4 };
        }
        v
    };
    let qseq0 = vec![qseq0_0, qseq0_1];

    // squeeze_a
    let mut order: Vec<(i32, usize)> = hits.iter().enumerate().map(|(i, h)| (h.a_off, i)).collect();
    order.sort_by_key(|&(a_off, _)| a_off);
    let mut as_off = 0i32;
    for &(_, hi) in &order {
        let r = &mut hits[hi];
        if r.a_off != as_off {
            for j in 0..r.cnt as usize {
                a[as_off as usize + j] = a[r.a_off as usize + j];
            }
            r.a_off = as_off;
        }
        as_off += r.cnt;
    }

    let n_regs = hits.len();
    let mut i = 0usize;
    while i < n_regs {
        let mut r2 = Hit::default();
        align1(opt, idx, qlen, &qseq0, &mut hits[i], &mut r2, a);
        if std::env::var("MB_DEBUG_ALIGN").is_ok() {
            if let Some(p) = &hits[i].p {
                let cig: String = p.cigar.iter()
                    .map(|c| format!("{}{}", c>>4, CIGAR_STR[(c&0xf) as usize] as char))
                    .collect();
                eprintln!("DBG hit{} qs={} qe={} ts={} te={} cigar={} cnt={}",
                    i, hits[i].qs, hits[i].qe, hits[i].ts, hits[i].te, cig, hits[i].cnt);
            }
        }
        if r2.cnt > 0 {
            hits.insert(i + 1, r2);
        }
        i += 1;
    }

    // filter_hits
    hits.retain(|r| {
        if r.p.is_none() { return r.cnt >= 0 && !r.flt; } // 保留无 p 的 hit（不应出现）
        let flt = r.flt || r.mlen < opt.min_chain_score || r.p.as_ref().unwrap().dp_max < opt.min_dp_max * opt.a;
        !flt
    });

    if !opt.is_sr_mode(qlen) {
        update_dp_max(qlen, hits, 0.9, opt.a, opt.b);
        hits.retain(|r| {
            if r.p.is_none() { return !r.flt; }
            let flt = r.flt || r.mlen < opt.min_chain_score || r.p.as_ref().unwrap().dp_max < opt.min_dp_max * opt.a;
            !flt
        });
    }

    // hit_sort: 按 dp_max（或 score）降序
    hits.sort_by(|a, b| {
        let sa = a.p.as_ref().map_or(a.score, |p| p.dp_max);
        let sb = b.p.as_ref().map_or(b.score, |p| p.dp_max);
        sb.cmp(&sa).then(a.hash.cmp(&b.hash))
    });
}

fn update_dp_max(qlen: i32, regs: &mut [Hit], frac: f64, a: i32, b: i32) {
    if regs.len() < 2 { return; }
    let mut max = -1i32; let mut max2 = -1i32;
    let mut max_i = -1isize; let mut max2_i = -1isize;
    for (i, r) in regs.iter().enumerate() {
        if let Some(p) = &r.p {
            if p.dp_max > max { max2 = max; max2_i = max_i; max = p.dp_max; max_i = i as isize; }
            else if p.dp_max > max2 { max2 = p.dp_max; max2_i = i as isize; }
        }
    }
    if max_i < 0 || max2_i < 0 { return; }
    let mi = max_i as usize;
    let m2i = max2_i as usize;
    if regs[mi].qe - regs[mi].qs < (qlen as f64 * frac) as i32 { return; }
    if regs[m2i].qe - regs[m2i].qs < ((regs[mi].qe - regs[mi].qs) as f64 * frac.sqrt()) as i32 { return; }
    let identity = event_identity(&regs[mi]);
    let mut div = 1.0 - identity;
    if div < 0.02 { div = 0.02; }
    let mut b2 = 0.5 / div;
    if (b2 * a as f64) < (b as f64) { b2 = a as f64 / b as f64; }
    for r in regs.iter_mut() {
        if r.p.is_none() { continue; }
        let new_dp = recal_max_dp(r, b2, a, qlen);
        r.p.as_mut().unwrap().dp_max = if new_dp < 0 { 0 } else { new_dp };
    }
}

fn event_identity(r: &Hit) -> f64 {
    if r.p.is_none() { return -1.0; }
    let (n_gap, n_gapo) = count_gaps(r);
    (r.mlen as f64) / (r.blen + r.p.as_ref().unwrap().n_ambi as i32 - n_gap + n_gapo) as f64
}

fn count_gaps(r: &Hit) -> (i32, i32) {
    if r.p.is_none() { return (-1, -1); }
    let mut n_gapo = 0; let mut n_gap = 0;
    for &c in &r.p.as_ref().unwrap().cigar {
        let op = c & 0xf;
        if op == CIGAR_INS || op == CIGAR_DEL { n_gapo += 1; n_gap += (c >> 4) as i32; }
    }
    (n_gap, n_gapo)
}

fn recal_max_dp(r: &Hit, b2: f64, match_sc: i32, qlen: i32) -> i32 {
    if r.p.is_none() { return -1; }
    let p = r.p.as_ref().unwrap();
    let mut n_gap = 0i32;
    let mut gap_cost = 0.0f64;
    for &c in &p.cigar {
        let op = c & 0xf; let len = (c >> 4) as i32;
        if op == CIGAR_INS || op == CIGAR_DEL {
            gap_cost += b2 * log2_approx(1.0 + len as f32) as f64;
            n_gap += len;
        }
    }
    let mut n_mis = r.blen + p.n_ambi as i32 - r.mlen - n_gap;
    n_mis += ((qlen - (r.qe - r.qs)) as f64 / b2 + 0.499) as i32;
    (match_sc as f64 * (r.mlen as f64 - b2 * n_mis as f64 - gap_cost) + 0.499) as i32
}
