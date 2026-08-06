//! 链化 DP，逐行对照 C 版 `lchain.c`（minimap2 风格，适配变长种子）。

use crate::l2bit::L2b;
use crate::seed::Anchor;

/// 两个锚点间的链化得分（对应 `comput_sc`）。无效返回 i32::MIN。
#[inline]
fn comput_sc(ai: &Anchor, aj: &Anchor, max_dist_x: i32, max_dist_y: i32, bw: i32, chn_pen_gap: f32) -> i32 {
    let dq = ai.qpos as i64 - aj.qpos as i64;
    if dq <= 0 || dq > max_dist_y as i64 + ai.len as i64 {
        return i32::MIN;
    }
    if ai.sid != aj.sid {
        return i32::MIN;
    }
    let dr = ai.tpos - aj.tpos;
    if dr <= 0 || dq > max_dist_x as i64 + ai.len as i64 {
        return i32::MIN;
    }
    let dd = if dr > dq { dr - dq } else { dq - dr };
    if dd > bw as i64 {
        return i32::MIN;
    }
    let dg = if dr < dq { dr } else { dq };
    let mut sc = if (ai.len as i64) < dg { ai.len as i64 } else { dg };
    if dd > 0 {
        let lin_pen = chn_pen_gap * dd as f32;
        let log_pen = crate::dna::log2_approx((dd + 1) as f32);
        sc -= (lin_pen + 0.5 * log_pen) as i64;
    }
    sc as i32
}

/// 回溯确定一条链的终点（对应 `mb_chain_bk_end`）。
fn chain_bk_end(max_drop: i32, z: &[(u64, i64)], f: &[i32], p: &[i64], t: &mut [i32], k: usize) -> i64 {
    let mut i = z[k].1;
    let mut end_i;
    let mut max_i = i;
    let mut max_s = 0i32;
    if i < 0 || t[i as usize] != 0 {
        return i;
    }
    loop {
        t[i as usize] = 2;
        end_i = p[i as usize];
        i = end_i;
        let s = if i < 0 {
            z[k].0 as i32
        } else {
            z[k].0 as i32 - f[i as usize]
        };
        if s > max_s {
            max_s = s;
            max_i = i;
        } else if max_s - s > max_drop {
            break;
        }
        if i < 0 || t[i as usize] != 0 {
            break;
        }
    }
    // 复位被标记的 t[]
    let mut i = z[k].1;
    while i >= 0 && i != end_i {
        t[i as usize] = 0;
        i = p[i as usize];
    }
    max_i
}

/// 回溯收集链（对应 `mb_chain_backtrack`）。
/// 返回 u（score<<32 | 锚点数）与 v（链内锚点下标，按链序拼接）。
#[allow(clippy::needless_range_loop)]
fn chain_backtrack(
    n: usize,
    f: &[i32],
    p: &[i64],
    v: &mut [i32],
    t: &mut [i32],
    min_sc: i32,
    max_drop: i32,
) -> Option<Vec<u64>> {
    let mut z: Vec<(u64, i64)> = Vec::new();
    for i in 0..n {
        if f[i] >= min_sc {
            z.push((f[i] as u64, i as i64));
        }
    }
    if z.is_empty() {
        return None;
    }
    // 基数排序等价：按 x 稳定升序
    z.sort_by_key(|&(x, _)| x);

    let mut u: Vec<u64> = Vec::new();
    t[..n].fill(0);
    let mut n_v = 0usize;
    for k in (0..z.len()).rev() {
        if t[z[k].1 as usize] == 0 {
            let n_v0 = n_v;
            let end_i = chain_bk_end(max_drop, &z, f, p, t, k);
            let mut i = z[k].1;
            while i != end_i {
                v[n_v] = i as i32;
                n_v += 1;
                t[i as usize] = 1;
                i = p[i as usize];
            }
            let sc = if i < 0 {
                z[k].0 as i32
            } else {
                z[k].0 as i32 - f[i as usize]
            };
            if sc >= min_sc && n_v > n_v0 {
                u.push((sc as u64) << 32 | (n_v - n_v0) as u64);
            } else {
                n_v = n_v0;
            }
        }
    }
    Some(u)
}

/// 链化主 DP（对应 `mb_lchain_dp`）。
///
/// 输入锚点须按拼接 tpos 有序。返回 (锚点数组重排后的链锚点, u)；
/// u[i] = score<<32 | 锚点数，链按首锚点目标位置有序（非得分序）。
#[allow(clippy::needless_range_loop)]
pub fn lchain_dp(
    l2b: &L2b,
    max_dist_x: i32,
    max_dist_y: i32,
    bw: i32,
    max_skip: i32,
    max_iter: i32,
    min_sc: i32,
    chn_pen_gap: f32,
    a: Vec<Anchor>,
) -> Option<(Vec<Anchor>, Vec<u64>)> {
    let n = a.len();
    if n == 0 {
        return None;
    }
    let max_dist_x = max_dist_x.max(bw);
    let max_dist_y = max_dist_y.max(bw);
    let max_drop = bw;

    let mut f = vec![0i32; n];
    let mut p = vec![0i64; n];
    let mut v = vec![0i32; n];
    let mut t = vec![0i32; n];
    let mut mmax_f = 0i32;
    let mut max_ii = -1i64;

    for i in 0..n {
        let mut max_j = -1i64;
        let mut max_f = a[i].len;
        let mut n_skip = 0i32;
        let mut j = i as i64 - 1;
        let j_floor = i as i64 - max_iter as i64;
        while j >= 0 && j >= j_floor {
            if a[i].tpos - a[j as usize].tpos >= max_dist_x as i64 + a[i].len as i64 {
                break;
            }
            let sc = comput_sc(&a[i], &a[j as usize], max_dist_x, max_dist_y, bw, chn_pen_gap);
            if sc != i32::MIN {
                let sc = sc + f[j as usize];
                if sc > max_f {
                    max_f = sc;
                    max_j = j;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                } else if t[j as usize] == i as i32 {
                    n_skip += 1;
                    if n_skip > max_skip {
                        break;
                    }
                }
            }
            if p[j as usize] >= 0 {
                t[p[j as usize] as usize] = i as i32;
            }
            j -= 1;
        }
        let end_j = j;
        if max_ii < 0 || a[i].tpos - a[max_ii as usize].tpos > max_dist_x as i64 + a[i].len as i64 {
            let mut max = i32::MIN;
            max_ii = -1;
            let mut jj = i as i64 - 1;
            while jj >= end_j && jj >= 0 {
                if max < f[jj as usize] {
                    max = f[jj as usize];
                    max_ii = jj;
                }
                jj -= 1;
            }
        }
        if max_ii >= 0 && max_ii < end_j {
            let tmp = comput_sc(&a[i], &a[max_ii as usize], max_dist_x, max_dist_y, bw, chn_pen_gap);
            if tmp != i32::MIN && max_f < tmp + f[max_ii as usize] {
                max_f = tmp + f[max_ii as usize];
                max_j = max_ii;
            }
        }
        f[i] = max_f;
        p[i] = max_j;
        v[i] = if max_j >= 0 && v[max_j as usize] > max_f {
            v[max_j as usize]
        } else {
            max_f
        };
        if max_ii < 0
            || (a[i].tpos - a[max_ii as usize].tpos <= max_dist_x as i64 + a[i].len as i64
                && f[max_ii as usize] < f[i])
        {
            max_ii = i as i64;
        }
        if mmax_f < max_f {
            mmax_f = max_f;
        }
    }

    let u = chain_backtrack(n, &f, &p, &mut v, &mut t, min_sc, max_drop)?;
    let n_u = u.len();
    if n_u == 0 {
        return None;
    }

    // compact_a：把每条链的锚点按序写入 b，再按目标位置排序链
    let n_v: usize = u.iter().map(|&x| (x & 0xffff_ffff) as usize).sum();
    let mut b: Vec<Anchor> = Vec::with_capacity(n_v);
    let mut k = 0usize;
    for &ui in &u {
        let ni = (ui & 0xffff_ffff) as usize;
        let k0 = k;
        for j in 0..ni {
            b.push(a[v[k0 + (ni - j - 1)] as usize]);
        }
        k += ni;
    }

    // 按首锚点拼接坐标排序链
    let mut w: Vec<(u64, u64)> = Vec::with_capacity(n_u);
    let mut k = 0usize;
    for (i, &ui) in u.iter().enumerate() {
        let first = &b[k];
        let ctg = &l2b.ctg[(first.sid >> 1) as usize];
        let x = (first.tpos + (ctg.off * 2 + ctg.len * (first.sid & 1) as u64) as i64) as u64;
        w.push((x, (k as u64) << 32 | i as u64));
        k += (ui & 0xffff_ffff) as usize;
    }
    w.sort_by_key(|&(x, _)| x); // 稳定，与 LSD 基数排序一致
    let mut u2 = vec![0u64; n_u];
    let mut a_out: Vec<Anchor> = Vec::with_capacity(n_v);
    for (i, &(_, wy)) in w.iter().enumerate() {
        let j = (wy & 0xffff_ffff) as usize;
        let off = (wy >> 32) as usize;
        let cnt = (u[j] & 0xffff_ffff) as usize;
        u2[i] = u[j];
        a_out.extend_from_slice(&b[off..off + cnt]);
    }
    Some((a_out, u2))
}
