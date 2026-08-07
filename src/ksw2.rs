//! ksw2 DP 比对的标量移植，对照 C 版 `ksw2_extd2_sse.c` / `ksw2.h`。
//!
//! 直接 i32 带状 DP，计算与 SSE 版等价的双仿射间隙递推：
//!   H(i,j) = min(max(H(i-1,j-1)+s, E1, F1, E2, F2), sc_mch)
//!   E1(i,j) = max(H(i,j-1)-q-e, E1(i,j-1)-e, -(q+e))   [DEL]
//!   F1(i,j) = max(H(i-1,j)-q-e, F1(i-1,j)-e, -(q+e))   [INS]
//!   E2/F2 同理用 q2/e2
//! 回溯矩阵 p 按 (query, target) 布局，backtrack 在 (i,j) 空间进行。

pub const NEG_INF: i32 = -0x40000000;

pub mod flag {
    pub const SCORE_ONLY: u32 = 0x01;
    pub const RIGHT: u32 = 0x02;
    pub const GENERIC_SC: u32 = 0x04;
    pub const APPROX_MAX: u32 = 0x08;
    pub const APPROX_DROP: u32 = 0x10;
    pub const EXTZ_ONLY: u32 = 0x40;
    pub const REV_CIGAR: u32 = 0x80;
}

pub const CIGAR_MATCH: u32 = 0;
pub const CIGAR_INS: u32 = 1;
pub const CIGAR_DEL: u32 = 2;

#[derive(Debug, Clone, Default)]
pub struct Extz {
    pub max: i32, pub zdropped: bool, pub max_q: i32, pub max_t: i32,
    pub mqe: i32, pub mqe_t: i32, pub mte: i32, pub mte_q: i32,
    pub score: i32, pub cigar: Vec<u32>, pub reach_end: bool,
}

pub fn reset_extz(ez: &mut Extz) {
    *ez = Extz { max: 0, score: NEG_INF, mqe: NEG_INF, mte: NEG_INF,
        max_q: -1, max_t: -1, mqe_t: -1, mte_q: -1, ..Default::default() };
}

fn push_cigar(c: &mut Vec<u32>, op: u32, len: u32) {
    if let Some(l) = c.last_mut() { if *l & 0xf == op { *l += len << 4; return; } }
    c.push(len << 4 | op);
}

pub fn gen_nt4_mat(mat: &mut [i8; 25], a: i8, b: i8, b_ts: i8, b_ambi: i8, mt: i32) {
    const M: usize = 5;
    let (a, b) = (a.abs(), if b > 0 { -b } else { b });
    let b_ambi = if b_ambi > 0 { -b_ambi } else { b_ambi };
    for i in 0..M - 1 {
        for j in 0..M - 1 { mat[i * M + j] = if i == j { a } else { b }; }
        mat[i * M + M - 1] = b_ambi;
    }
    for j in 0..M { mat[(M - 1) * M + j] = 0; }
    if b_ts != 0 && mt == 0 { mat[M + 3] = -b_ts; mat[3 * M + 1] = -b_ts; }
}

fn backtrack(is_rev: bool, p: &[Vec<u8>], i0: i32, j0: i32, cigar: &mut Vec<u32>) {
    let (mut i, mut j, mut state) = (i0, j0, 0u8);
    while i >= 0 && j >= 0 {
        let tmp = p[i as usize][j as usize];
        if state == 0 { state = tmp & 7; }
        else if tmp >> (state + 2) & 1 == 0 { state = 0; }
        if state == 0 { state = tmp & 7; }
        match state {
            0 => { push_cigar(cigar, CIGAR_MATCH, 1); i -= 1; j -= 1; }
            1 | 3 => { push_cigar(cigar, CIGAR_DEL, 1); j -= 1; }
            _ => { push_cigar(cigar, CIGAR_INS, 1); i -= 1; }
        }
    }
    if i >= 0 { push_cigar(cigar, CIGAR_INS, (i + 1) as u32); }
    if j >= 0 { push_cigar(cigar, CIGAR_DEL, (j + 1) as u32); }
    if !is_rev { cigar.reverse(); }
}

#[allow(clippy::too_many_arguments)]
pub fn ksw_extd2(
    qlen: i32, query: &[u8], tlen: i32, target: &[u8], m: i32, mat: &[i8],
    mut q: i8, mut e: i8, mut q2: i8, mut e2: i8, mut w: i32,
    zdrop: i32, end_bonus: i32, flags: u32, ez: &mut Extz,
) {
    reset_extz(ez);
    if m <= 1 || qlen <= 0 || tlen <= 0 { return; }
    if (q2 as i32 + e2 as i32) < (q as i32 + e as i32) {
        std::mem::swap(&mut q, &mut q2); std::mem::swap(&mut e, &mut e2);
    }
    let with_cigar = flags & flag::SCORE_ONLY == 0;
    let right = flags & flag::RIGHT != 0;
    let generic = flags & flag::GENERIC_SC != 0;
    let qlen = qlen as usize; let tlen = tlen as usize; let m = m as usize;
    let (qe, qe2) = (q as i32 + e as i32, q2 as i32 + e2 as i32);
    let sc_mch = mat[0] as i32;
    let sc_n = if mat[m * m - 1] == 0 { -e2 as i32 } else { mat[m * m - 1] as i32 };
    let m1 = (m - 1) as u8;

    if w < 0 { w = tlen.max(qlen) as i32; }
    let w = w as usize;
    let mut max_sc = mat[0] as i32; let mut min_sc = mat[1] as i32;
    for &v in mat.iter().take(m * m).skip(1) {
        if (v as i32) > max_sc { max_sc = v as i32; }
        if (v as i32) < min_sc { min_sc = v as i32; }
    }
    if -min_sc > 2 * qe { return; }
    if qlen as i64 * tlen as i64 > 100_000_000 { ez.zdropped = true; return; }

    // 行滚动数组
    let mut hp = vec![NEG_INF; tlen];  // H(i-1, j)
    let mut hc = vec![NEG_INF; tlen];  // H(i, j)
    let mut e1 = vec![NEG_INF; tlen];  // E1(i, j) — 同行左→右
    let mut e2arr = vec![NEG_INF; tlen];
    let mut f1p = vec![NEG_INF; tlen]; // F1(i-1, j)
    let mut f1c = vec![NEG_INF; tlen]; // F1(i, j)
    let mut f2p = vec![NEG_INF; tlen];
    let mut f2c = vec![NEG_INF; tlen];

    let mut p: Vec<Vec<u8>> = if with_cigar { vec![vec![0u8; tlen]; qlen] } else { Vec::new() };

    for i in 0..qlen {
        let j_lo = i.saturating_sub(w);
        let j_hi = (i + w).min(tlen - 1);
        for j in 0..tlen { hc[j] = NEG_INF; e1[j] = NEG_INF; e2arr[j] = NEG_INF; f1c[j] = NEG_INF; f2c[j] = NEG_INF; }

        for j in j_lo..=j_hi {
            let s = if generic {
                mat[target[j] as usize * m + query[i] as usize] as i32
            } else {
                let (sq, qt) = (target[j], query[i]);
                if sq == m1 || qt == m1 { sc_n } else if sq == qt { sc_mch } else { mat[1] as i32 }
            };

            // 对角线 H(i-1, j-1)
            let h_diag = if i == 0 && j == 0 { 0 } else if j > 0 { hp[j - 1] } else { NEG_INF };
            let mut h = h_diag.saturating_add(s);
            let mut d = 0u8;

            // E1 (DEL): H(i, j-1) 开 / E1(i, j-1) 延伸 / 地板 -(q+e)
            if j > 0 {
                let open = hc[j - 1].saturating_sub(qe);
                let ext = e1[j - 1].saturating_sub(e as i32);
                let val = open.max(ext).max(-qe);
                e1[j] = val;
                // 决策：left 用 a > h（严格），right 用 a >= h
                let a_raw = (hc[j - 1].saturating_sub(q as i32)).max(e1[j - 1]); // 地板前候选
                if !right {
                    if val > h { d = 1; h = val; }
                } else if val >= h { d = 1; h = val; }
                if !right { if a_raw > 0 { d |= 0x08; } }
                else { if a_raw >= 0 { d |= 0x08; } }
            }

            // F1 (INS): H(i-1, j) 开 / F1(i-1, j) 延伸 / 地板
            {
                let open = hp[j].saturating_sub(qe);
                let ext = f1p[j].saturating_sub(e as i32);
                let val = open.max(ext).max(-qe);
                f1c[j] = val;
                let a_raw = (hp[j].saturating_sub(q as i32)).max(f1p[j]);
                if !right {
                    if val > h { d = 2; h = val; }
                } else if val >= h { d = 2; h = val; }
                if !right { if a_raw > 0 { d |= 0x10; } }
                else { if a_raw >= 0 { d |= 0x10; } }
            }

            // E2 (long DEL)
            if j > 0 {
                let open = hc[j - 1].saturating_sub(qe2);
                let ext = e2arr[j - 1].saturating_sub(e2 as i32);
                let val = open.max(ext).max(-qe2);
                e2arr[j] = val;
                let a_raw = (hc[j - 1].saturating_sub(q2 as i32)).max(e2arr[j - 1]);
                if !right {
                    if val > h { d = 3; h = val; }
                } else if val >= h { d = 3; h = val; }
                if !right { if a_raw > 0 { d |= 0x20; } }
                else { if a_raw >= 0 { d |= 0x20; } }
            }

            // F2 (long INS)
            {
                let open = hp[j].saturating_sub(qe2);
                let ext = f2p[j].saturating_sub(e2 as i32);
                let val = open.max(ext).max(-qe2);
                f2c[j] = val;
                let a_raw = (hp[j].saturating_sub(q2 as i32)).max(f2p[j]);
                if !right {
                    if val > h { d = 4; h = val; }
                } else if val >= h { d = 4; h = val; }
                if !right { if a_raw > 0 { d |= 0x40; } }
                else { if a_raw >= 0 { d |= 0x40; } }
            }

            hc[j] = h;
            if with_cigar { p[i][j] = d; }

            // max / zdrop：APPROX_MAX && !APPROX_DROP 时跳过（C 版同样不更新 ez.max）
            let track_max = !(flags & flag::APPROX_MAX != 0 && flags & flag::APPROX_DROP == 0);
            if track_max {
                let (r, t) = ((i + j) as i32, j as i32);
                if h > ez.max {
                    ez.max = h; ez.max_t = t; ez.max_q = r - t;
                } else if t >= ez.max_t && r - t >= ez.max_q {
                    let (tl, ql) = (t - ez.max_t, (r - t) - ez.max_q);
                    let l = if tl > ql { tl - ql } else { ql - tl };
                    if zdrop >= 0 && ez.max - h > zdrop + l * e2 as i32 { ez.zdropped = true; }
                }
                if j == tlen - 1 && h > ez.mte { ez.mte = h; ez.mte_q = i as i32; }
                if i == qlen - 1 && h > ez.mqe { ez.mqe = h; ez.mqe_t = j as i32; }
            }
        }

        std::mem::swap(&mut hp, &mut hc);
        std::mem::swap(&mut f1p, &mut f1c);
        std::mem::swap(&mut f2p, &mut f2c);
        if ez.zdropped { break; }
    }

    if !ez.zdropped { ez.score = hp[tlen - 1]; }

    if with_cigar {
        let rev = flags & flag::REV_CIGAR != 0;
        if !ez.zdropped && flags & flag::EXTZ_ONLY == 0 {
            backtrack(rev, &p, qlen as i32 - 1, tlen as i32 - 1, &mut ez.cigar);
        } else if flags & flag::EXTZ_ONLY != 0 && ez.mqe + end_bonus > ez.max {
            ez.reach_end = true;
            backtrack(rev, &p, qlen as i32 - 1, ez.mqe_t, &mut ez.cigar);
        } else if ez.max_q >= 0 && ez.max_t >= 0 {
            backtrack(rev, &p, ez.max_q - 1, ez.max_t, &mut ez.cigar);
        }
    }
}

pub fn ll_local(query: &[u8], target: &[u8], m: usize, mat: &[i8], gapo: i32, gape: i32) -> i32 {
    let (qlen, tlen) = (query.len(), target.len());
    if qlen == 0 || tlen == 0 { return 0; }
    let mut h = vec![0i32; qlen + 1];
    let mut er = vec![NEG_INF; qlen + 1];
    let mut best = 0i32;
    for j in 1..=tlen {
        let tj = target[j - 1] as usize;
        let mut f = NEG_INF;
        let mut hd = 0i32;
        for i in 1..=qlen {
            let s = mat[tj * m + query[i - 1] as usize] as i32;
            let hp = h[i];
            er[i] = (er[i] - gape).max(hd - gapo - gape);
            f = (f - gape).max(hp - gapo - gape);
            let mut hh = hd + s;
            if er[i] > hh { hh = er[i]; }
            if f > hh { hh = f; }
            if hh < 0 { hh = 0; }
            h[i] = hh;
            if hh > best { best = hh; }
            hd = hp;
        }
    }
    best
}
