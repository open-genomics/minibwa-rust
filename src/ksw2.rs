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
    pub max: i32,
    pub zdropped: bool,
    pub max_q: i32,
    pub max_t: i32,
    pub mqe: i32,
    pub mqe_t: i32,
    pub mte: i32,
    pub mte_q: i32,
    pub score: i32,
    pub cigar: Vec<u32>,
    pub reach_end: bool,
}

pub fn reset_extz(ez: &mut Extz) {
    *ez = Extz {
        max: 0,
        score: NEG_INF,
        mqe: NEG_INF,
        mte: NEG_INF,
        max_q: -1,
        max_t: -1,
        mqe_t: -1,
        mte_q: -1,
        ..Default::default()
    };
}

pub fn push_cigar_vec(c: &mut Vec<u32>, op: u32, len: u32) {
    if let Some(l) = c.last_mut() {
        if *l & 0xf == op {
            *l += len << 4;
            return;
        }
    }
    c.push(len << 4 | op);
}

pub fn gen_nt4_mat(mat: &mut [i8; 25], a: i8, b: i8, b_ts: i8, b_ambi: i8, mt: i32) {
    const M: usize = 5;
    let (a, b) = (a.abs(), if b > 0 { -b } else { b });
    let b_ambi = if b_ambi > 0 { -b_ambi } else { b_ambi };
    for i in 0..M - 1 {
        for j in 0..M - 1 {
            mat[i * M + j] = if i == j { a } else { b };
        }
        mat[i * M + M - 1] = b_ambi;
    }
    for j in 0..M {
        mat[(M - 1) * M + j] = b_ambi;
    }
    if b_ts != 0 && b_ts != b {
        let b_ts = if b_ts > 0 { -b_ts } else { b_ts };
        mat[2] = b_ts;
        mat[M + 3] = b_ts;
        mat[2 * M] = b_ts;
        mat[3 * M + 1] = b_ts;
    }
    if mt == 1 {
        mat[M + 3] = a;
    } else if mt == 2 {
        mat[2 * M] = a;
    }
}

#[allow(clippy::too_many_arguments)]
fn ksw_extd2_antidiag(
    qlen: usize,
    query: &[u8],
    tlen: usize,
    target: &[u8],
    m: usize,
    mat: &[i8],
    q: i8,
    e: i8,
    q2: i8,
    e2: i8,
    w: i32,
    zdrop: i32,
    end_bonus: i32,
    flags: u32,
    qe_for_h: i32,
    with_cigar: bool,
    right: bool,
    ez: &mut Extz,
) {
    let generic = flags & flag::GENERIC_SC != 0;
    let approx_max = flags & flag::APPROX_MAX != 0;
    let extz_only = flags & flag::EXTZ_ONLY != 0;
    let approx_drop = flags & flag::APPROX_DROP != 0;
    let qe = q as i32 + e as i32;
    let qe2 = q2 as i32 + e2 as i32;
    let sc_mch = mat[0];
    let sc_n = if mat[m * m - 1] == 0 { -e2 } else { mat[m * m - 1] };
    let m1 = (m - 1) as u8;
    let max_sc = mat.iter().take(m * m).map(|&v| v as i32).max().unwrap();
    let tpad = tlen.div_ceil(16) * 16;

    let mut u = vec![-(qe as i8); tpad];
    let mut v = vec![-(qe as i8); tpad];
    let mut x = vec![-(qe as i8); tpad];
    let mut y = vec![-(qe as i8); tpad];
    let mut x2 = vec![-(qe2 as i8); tpad];
    let mut y2 = vec![-(qe2 as i8); tpad];
    let mut scores = vec![0i8; tpad];
    let mut h = if approx_max { None } else { Some(vec![NEG_INF; tlen]) };
    let mut h0 = 0i32;
    let mut last_h0_t = 0i32;
    let mut last_max_h = 0i32;
    let mut last_st = -1i32;
    let mut last_en = -1i32;
    // 回溯矩阵按反斜角存成扁平缓冲：每行一个 `off`/`off_end` 区间，行首偏移记在
    // `p_off` 里。逐行 `Vec` 会为每个反斜角各分配一次（长读动辄上千次小分配），
    // 扁平缓冲只在容量不足时扩容。
    let mut p: Vec<u8> = Vec::new();
    let mut p_off: Vec<usize> = Vec::new();
    let mut off = Vec::new();
    let mut off_end = Vec::new();
    let long_thres = if e != e2 {
        (q2 as i32 - q as i32) / (e as i32 - e2 as i32) - 1
    } else {
        0
    };
    let long_thres = if q2 as i32 + e2 as i32 + long_thres * e2 as i32 > q as i32 + e as i32 + long_thres * e as i32 {
        long_thres + 1
    } else {
        long_thres
    };
    let long_diff = long_thres * (e as i32 - e2 as i32) - (q2 as i32 - q as i32) - e2 as i32;
    let band = if w < 0 { qlen.max(tlen) as i32 } else { w };

    for r in 0..(qlen + tlen - 1) as i32 {
        let mut st0 = 0i32;
        let mut en0 = tlen as i32 - 1;
        st0 = st0.max(r - qlen as i32 + 1);
        en0 = en0.min(r);
        st0 = st0.max((r - band + 1) >> 1);
        en0 = en0.min((r + band) >> 1);
        if st0 > en0 {
            ez.zdropped = true;
            break;
        }

        let st = (st0 / 16) * 16;
        let en = ((en0 + 16) / 16) * 16 - 1;
        let mut drow = if with_cigar {
            vec![0u8; (en - st + 1) as usize]
        } else {
            Vec::new()
        };

        let mut x1 = -(qe as i8);
        let mut x21 = -(qe2 as i8);
        let mut v1 = -(qe as i8);
        if st > 0 {
            if st > last_st && st - 1 <= last_en {
                let t = (st - 1) as usize;
                x1 = x[t];
                x21 = x2[t];
                v1 = v[t];
            }
        } else {
            v1 = if r == 0 {
                -(qe as i8)
            } else if r < long_thres {
                -e
            } else if r == long_thres {
                long_diff as i8
            } else {
                -e2
            };
        }
        if en >= r {
            let t = r as usize;
            y[t] = -(qe as i8);
            y2[t] = -(qe2 as i8);
            u[t] = if r == 0 {
                -(qe as i8)
            } else if r < long_thres {
                -e
            } else if r == long_thres {
                long_diff as i8
            } else {
                -e2
            };
        }

        if generic {
            for t_i in st0..=en0 {
                let t = t_i as usize;
                scores[t] = mat[target[t] as usize * m + query[(r - t_i) as usize] as usize];
            }
        } else {
            let mut start = st0;
            while start <= en0 {
                for lane in 0..16i32 {
                    let t_i = start + lane;
                    let t = t_i as usize;
                    if t >= tpad {
                        continue;
                    }
                    let target_base = if t < tlen { target[t] } else { 0 };
                    let qrr_pos = qlen as i32 - 1 - r + t_i;
                    let query_base = if qrr_pos >= 0 && qrr_pos < qlen as i32 {
                        query[qlen - 1 - qrr_pos as usize]
                    } else if qrr_pos < 0 {
                        let sf_pos = tpad as i32 + qrr_pos;
                        if sf_pos >= 0 && sf_pos < tpad as i32 && (sf_pos as usize) < tlen {
                            target[sf_pos as usize]
                        } else {
                            0
                        }
                    } else {
                        0
                    };
                    scores[t] = if target_base == m1 || query_base == m1 {
                        sc_n
                    } else if target_base == query_base {
                        sc_mch
                    } else {
                        mat[1]
                    };
                }
                start += 16;
            }
        }

        for t_i in st..=en {
            let t = t_i as usize;
            let old_x = x[t];
            let old_x2 = x2[t];
            let old_v = v[t];
            let old_y = y[t];
            let old_y2 = y2[t];
            let old_u = u[t];
            let score = scores[t];

            let a = x1.wrapping_add(v1);
            let b = old_y.wrapping_add(old_u);
            let a2 = x21.wrapping_add(v1);
            let b2 = old_y2.wrapping_add(old_u);
            let mut z = score;
            let mut d = 0u8;
            if !right {
                if a > z {
                    d = 1;
                    z = a;
                }
                if b > z {
                    d = 2;
                    z = b;
                }
                if a2 > z {
                    d = 3;
                    z = a2;
                }
                if b2 > z {
                    d = 4;
                    z = b2;
                }
            } else {
                if z <= a {
                    d = 1;
                    z = a;
                }
                if z <= b {
                    d = 2;
                    z = b;
                }
                if z <= a2 {
                    d = 3;
                    z = a2;
                }
                if z <= b2 {
                    d = 4;
                    z = b2;
                }
            }
            z = z.min(sc_mch);

            u[t] = z.wrapping_sub(v1);
            v[t] = z.wrapping_sub(old_u);

            let z_q = z.wrapping_sub(q);
            let a = a.wrapping_sub(z_q);
            let b = b.wrapping_sub(z_q);
            let z_q2 = z.wrapping_sub(q2);
            let a2 = a2.wrapping_sub(z_q2);
            let b2 = b2.wrapping_sub(z_q2);
            let positive = |v: i8| if right { v >= 0 } else { v > 0 };
            x[t] = if positive(a) { a } else { 0 }.wrapping_sub(qe as i8);
            y[t] = if positive(b) { b } else { 0 }.wrapping_sub(qe as i8);
            x2[t] = if positive(a2) { a2 } else { 0 }.wrapping_sub(qe2 as i8);
            y2[t] = if positive(b2) { b2 } else { 0 }.wrapping_sub(qe2 as i8);
            if with_cigar {
                if positive(a) {
                    d |= 0x08;
                }
                if positive(b) {
                    d |= 0x10;
                }
                if positive(a2) {
                    d |= 0x20;
                }
                if positive(b2) {
                    d |= 0x40;
                }
                drow[(t_i - st) as usize] = d;
            }

            x1 = old_x;
            x21 = old_x2;
            v1 = old_v;
        }

        if with_cigar {
            off.push(st);
            off_end.push(en);
            p_off.push(p.len());
            p.extend_from_slice(&drow);
        }

        if approx_max {
            if r > 0 {
                if last_h0_t >= st0 && last_h0_t <= en0 && last_h0_t + 1 >= st0 && last_h0_t < en0 {
                    let d0 = v[last_h0_t as usize] as i32;
                    let d1 = u[(last_h0_t + 1) as usize] as i32;
                    if d0 > d1 {
                        h0 += d0;
                    } else {
                        h0 += d1;
                        last_h0_t += 1;
                    }
                } else if last_h0_t >= st0 && last_h0_t <= en0 {
                    h0 += v[last_h0_t as usize] as i32;
                } else {
                    last_h0_t += 1;
                    h0 += u[last_h0_t as usize] as i32;
                }
                if approx_drop && zdrop >= 0 {
                    let t = last_h0_t;
                    if h0 > ez.max {
                        ez.max = h0;
                        ez.max_t = t;
                        ez.max_q = r - t;
                    } else if t >= ez.max_t && r - t >= ez.max_q {
                        let l = (t - ez.max_t - (r - t - ez.max_q)).abs();
                        if ez.max - h0 > zdrop + l * e2 as i32 {
                            ez.zdropped = true;
                            break;
                        }
                    }
                }
                if r == (qlen + tlen - 2) as i32 && en0 == tlen as i32 - 1 {
                    ez.score = h0;
                }
            } else {
                h0 = v[0] as i32 - qe_for_h;
                last_h0_t = 0;
            }
        } else {
            let h = h.as_mut().unwrap();
            let (max_h, max_t) = if r == 0 {
                h[0] = v[0] as i32 - qe_for_h;
                (h[0], 0)
            } else {
                let en = en0 as usize;
                h[en] = if en0 > 0 {
                    h[en - 1] + u[en] as i32
                } else {
                    h[en] + v[en] as i32
                };
                let mut max_h = h[en];
                let mut max_t = en0;
                let mut lane_h = [max_h; 4];
                let mut lane_t = [max_t; 4];
                let en1 = st0 + (en0 - st0) / 4 * 4;
                let mut t_i = st0;
                while t_i < en1 {
                    for lane in 0usize..4 {
                        let t = t_i as usize + lane;
                        h[t] += v[t] as i32;
                        if h[t] > lane_h[lane] {
                            lane_h[lane] = h[t];
                            lane_t[lane] = t_i;
                        }
                    }
                    t_i += 4;
                }
                for lane in 0..4 {
                    if max_h < lane_h[lane] {
                        max_h = lane_h[lane];
                        max_t = lane_t[lane] + lane as i32;
                    }
                }
                while t_i < en0 {
                    let t = t_i as usize;
                    h[t] += v[t] as i32;
                    if h[t] > max_h {
                        max_h = h[t];
                        max_t = t_i;
                    }
                    t_i += 1;
                }
                (max_h, max_t)
            };
            if en0 == tlen as i32 - 1 && h[en0 as usize] > ez.mte {
                ez.mte = h[en0 as usize];
                ez.mte_q = r - en0;
            }
            if r - st0 == qlen as i32 - 1 && h[st0 as usize] > ez.mqe {
                ez.mqe = h[st0 as usize];
                ez.mqe_t = st0;
            }
            if max_h > ez.max {
                ez.max = max_h;
                ez.max_t = max_t;
                ez.max_q = r - max_t;
            } else if max_t >= ez.max_t && r - max_t >= ez.max_q {
                let l = (max_t - ez.max_t - (r - max_t - ez.max_q)).abs();
                if zdrop >= 0 && ez.max - max_h > zdrop + l * e2 as i32 {
                    ez.zdropped = true;
                    break;
                }
            }
            if r == (qlen + tlen - 2) as i32 && en0 == tlen as i32 - 1 {
                ez.score = h[tlen - 1];
            }
            if extz_only {
                let r_h = last_max_h.max(max_h);
                let rq = qlen as i32 - (r - st0);
                let rt = tlen as i32 - en0;
                let rm = if rq >= tlen as i32 - st0 {
                    tlen as i32 - st0
                } else if rt >= qlen as i32 - (r - en0) {
                    qlen as i32 - (r - en0)
                } else {
                    tlen as i32 + qlen as i32 - 1 - r
                };
                if r_h + rm * max_sc + end_bonus < ez.max {
                    break;
                }
            }
            last_max_h = max_h;
        }
        last_st = st;
        last_en = en;
    }

    if with_cigar {
        if !ez.zdropped && !extz_only {
            backtrack(
                flags & flag::REV_CIGAR != 0,
                &p,
                &p_off,
                &off,
                &off_end,
                tlen as i32 - 1,
                qlen as i32 - 1,
                &mut ez.cigar,
            );
        } else if extz_only && ez.mqe + end_bonus > ez.max {
            ez.reach_end = true;
            backtrack(
                flags & flag::REV_CIGAR != 0,
                &p,
                &p_off,
                &off,
                &off_end,
                ez.mqe_t,
                qlen as i32 - 1,
                &mut ez.cigar,
            );
        } else if ez.max_t >= 0 && ez.max_q >= 0 {
            backtrack(
                flags & flag::REV_CIGAR != 0,
                &p,
                &p_off,
                &off,
                &off_end,
                ez.max_t,
                ez.max_q,
                &mut ez.cigar,
            );
        }
    }
}

fn backtrack(
    is_rev: bool,
    p: &[u8],
    p_off: &[usize],
    off: &[i32],
    off_end: &[i32],
    i0: i32,
    j0: i32,
    cigar: &mut Vec<u32>,
) {
    let (mut i, mut j, mut state) = (i0, j0, 0u8);
    while i >= 0 && j >= 0 {
        let r = (i + j) as usize;
        let mut force_state = -1i32;
        if i < off[r] {
            force_state = 2;
        }
        if i > off_end[r] {
            force_state = 1;
        }
        let tmp = if force_state < 0 {
            p[p_off[r] + (i - off[r]) as usize]
        } else {
            0
        };
        if state == 0 {
            state = tmp & 7;
        } else if tmp >> (state + 2) & 1 == 0 {
            state = 0;
        }
        if state == 0 {
            state = tmp & 7;
        }
        if force_state >= 0 {
            state = force_state as u8;
        }
        match state {
            0 => {
                push_cigar_vec(cigar, CIGAR_MATCH, 1);
                i -= 1;
                j -= 1;
            }
            1 | 3 => {
                push_cigar_vec(cigar, CIGAR_DEL, 1);
                i -= 1;
            }
            _ => {
                push_cigar_vec(cigar, CIGAR_INS, 1);
                j -= 1;
            }
        }
    }
    if i >= 0 {
        push_cigar_vec(cigar, CIGAR_DEL, (i + 1) as u32);
    }
    if j >= 0 {
        push_cigar_vec(cigar, CIGAR_INS, (j + 1) as u32);
    }
    if !is_rev {
        cigar.reverse();
    }
}

#[allow(clippy::too_many_arguments)]
pub fn ksw_extd2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    mut q: i8,
    mut e: i8,
    mut q2: i8,
    mut e2: i8,
    mut w: i32,
    zdrop: i32,
    end_bonus: i32,
    flags: u32,
    ez: &mut Extz,
) {
    reset_extz(ez);
    if m <= 1 || qlen <= 0 || tlen <= 0 {
        return;
    }
    let qe_for_h = q as i32 + e as i32;
    if (q2 as i32 + e2 as i32) < (q as i32 + e as i32) {
        std::mem::swap(&mut q, &mut q2);
        std::mem::swap(&mut e, &mut e2);
    }
    let with_cigar = flags & flag::SCORE_ONLY == 0;
    let right = flags & flag::RIGHT != 0;
    let qlen = qlen as usize;
    let tlen = tlen as usize;
    let m = m as usize;
    let qe = q as i32 + e as i32;

    if w < 0 {
        w = tlen.max(qlen) as i32;
    }
    let mut min_sc = mat[1] as i32;
    for &v in mat.iter().take(m * m).skip(1) {
        if (v as i32) < min_sc {
            min_sc = v as i32;
        }
    }
    if -min_sc > 2 * qe {
        return;
    }
    if qlen as i64 * tlen as i64 > 100_000_000 {
        ez.zdropped = true;
        return;
    }

    ksw_extd2_antidiag(
        qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flags, qe_for_h, with_cigar, right, ez,
    );
}

/// 线性空间局部比对结果（对应 `ksw_llrst_t`）。
#[derive(Debug, Clone, Copy, Default)]
pub struct LlRst {
    pub score: i32,
    pub te: i32,
    pub qe: i32,
    pub score2: i32,
    pub te2: i32,
}

/// 线性空间局部比对打分（对应 `ksw_ll_u8_core` / `ksw_ll_i16_core` 的标量版）。
///
/// 语义与 C 版一致：维护 H/E/F 三状态，记录每行最大值与达到全局最大时的 H 行
/// （Hmax），据此给出 te（target 端）与 qe（query 端）；`minsc` 对应 xtra 低 16 位
/// （`KSW_LL_SUBO | opt->min_len`），用于记录高分区段并计算次优分 score2。
/// 采用 i32 无截断计算：缩放后 match=+1，u8 版在 max_sc < 251 时不会饱和，
/// i16 版无截断，故本实现与 C 版两条路径等价。
#[allow(clippy::needless_range_loop)]
pub fn ll_core(query: &[u8], target: &[u8], m: usize, mat: &[i8], gapo: i32, gape: i32, minsc: i32) -> LlRst {
    let (qlen, tlen) = (query.len(), target.len());
    let mut r = LlRst {
        score: 0,
        te: -1,
        qe: -1,
        score2: -1,
        te2: -1,
    };
    if qlen == 0 || tlen == 0 {
        return r;
    }
    let mut mat_max = 0i32;
    for &v in mat.iter().take(m * m) {
        if (v as i32) > mat_max {
            mat_max = v as i32;
        }
    }
    let mut h = vec![0i32; qlen + 1]; // H(i-1, ·)；h[j+1] 对应 query 位置 j
    let mut e = vec![0i32; qlen + 1]; // E(i, ·)（插入状态）
    let mut hmax: Vec<i32> = vec![0; qlen]; // 全局最大行 H(i, ·)，索引即 query 位置
    let mut gmax = 0i32;
    let mut te = -1i32;
    let mut b: Vec<u64> = Vec::new(); // 高 32 位：行最大分；低 32 位：行号
    for i in 0..tlen {
        let tj = target[i] as usize;
        let mut f = 0i32; // F(i, j)（删除状态），每行重置
        let mut imax = 0i32;
        let mut hd = 0i32; // H(i-1, j-1)
        for j in 0..qlen {
            let hp = h[j + 1]; // H(i-1, j)
            let s = mat[tj * m + query[j] as usize] as i32;
            let mut hh = hd + s;
            if e[j + 1] > hh {
                hh = e[j + 1];
            }
            if f > hh {
                hh = f;
            }
            if hh < 0 {
                hh = 0;
            }
            h[j + 1] = hh;
            e[j + 1] = (e[j + 1] - gape).max(hh - gapo - gape);
            f = (f - gape).max(hh - gapo - gape);
            if hh > imax {
                imax = hh;
            }
            hd = hp;
        }
        if imax >= minsc {
            if b.is_empty() || (b.last().unwrap() & 0xffff_ffff) as i64 + 1 != i as i64 {
                b.push((imax as u64) << 32 | i as u64);
            } else if (b.last().unwrap() >> 32) < imax as u64 {
                *b.last_mut().unwrap() = (imax as u64) << 32 | i as u64;
            }
        }
        if imax > gmax {
            gmax = imax;
            te = i as i32;
            hmax.copy_from_slice(&h[1..]);
        }
    }
    r.score = gmax;
    r.te = te;
    // qe：gmax 行的 H 最大值所在 query 位置（对应 C 遍历 Hmax）
    let mut max = -1i32;
    for (i, &v) in hmax.iter().enumerate() {
        let ii = i as i32;
        if v > max {
            max = v;
            r.qe = ii;
        } else if v == max && ii < r.qe {
            r.qe = ii;
        }
    }
    // score2：离 te 窗口外最大的区段分（窗口宽度用 max(mat)，对应 C 版 q->max）
    if !b.is_empty() {
        let w = if mat_max > 0 { (gmax + mat_max - 1) / mat_max } else { 0 };
        let low = te - w;
        let high = te + w;
        for &bv in &b {
            let e2 = (bv & 0xffff_ffff) as i32;
            if (e2 < low || e2 > high) && (bv >> 32) as i32 > r.score2 {
                r.score2 = (bv >> 32) as i32;
                r.te2 = e2;
            }
        }
    }
    r
}

/// 纯分数版线性空间局部比对（只取最优分）。
pub fn ll_local(query: &[u8], target: &[u8], m: usize, mat: &[i8], gapo: i32, gape: i32) -> i32 {
    // xtra=0 → minsc=0x10000，b[] 永不记录，仅返回最优分
    ll_core(query, target, m, mat, gapo, gape, 0x10000).score
}
