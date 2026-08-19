//! 索引构建：FASTA → `.l2b` + `.mbw`，逐行对照 C 版 `index.c` 的 libsais 路径。

use anyhow::{bail, Result};

use crate::bwt::Bwt;
use crate::l2bit::{open_fasta, L2b};

/// libsais 要求的 SA 数组前置空闲空间（与 C 版 `fs` 一致）。
const FS: usize = 10000;

/// 为参考基因组 FASTA 建索引（对应 `main_index` 的默认 libsais 路径）。
/// 生成 `<prefix>.l2b` 与 `<prefix>.mbw`；`is_meth` 时再写 `<prefix>.meth.mbw`。
pub fn build(fasta: &str, prefix: &str, seed: u64, sa_bit: u32, threads: usize, is_meth: bool) -> Result<()> {
    let reader = open_fasta(fasta)?;
    let l2b = L2b::import(reader, seed)?;
    if l2b.tot_len == 0 {
        bail!("参考基因组为空: {fasta}");
    }
    eprintln!("[M::index] {} 条序列，共 {} bp", l2b.ctg.len(), l2b.tot_len);
    let bwt = bwt_libsais(&l2b, sa_bit, true, false, threads)?;
    let fn_l2b = format!("{prefix}.l2b");
    let fn_mbw = format!("{prefix}.mbw");
    l2b.save(&fn_l2b)?;
    bwt.save(&fn_mbw)?;
    if is_meth {
        let bwt_meth = bwt_libsais(&l2b, sa_bit, true, true, threads)?;
        let fn_meth = format!("{prefix}.meth.mbw");
        bwt_meth.save(&fn_meth)?;
        eprintln!("[M::index] 索引已写入 {fn_l2b} 与 {fn_mbw} 与 {fn_meth}");
    } else {
        eprintln!("[M::index] 索引已写入 {fn_l2b} 与 {fn_mbw}");
    }
    Ok(())
}

/// libsais 构建 BWT + 采样 SA（对应 `mb_bwt_libsais`）。
pub(crate) fn bwt_libsais(l2b: &L2b, sa_bit: u32, both_strand: bool, is_meth: bool, threads: usize) -> Result<Bwt> {
    let copies = if is_meth { 2u64 } else { 1 };
    let strands = if both_strand { 2u64 } else { 1 };
    let len = l2b.tot_len * copies * strands;
    let use_int32 = len + (FS as u64) < (i32::MAX as u64);

    // 拼接文本：对照 `mb_bwt_libsais` 的 is_meth / both_strand 分支
    let mut seq = vec![0u8; len as usize];
    let mut j = 0usize;
    if is_meth {
        for i in 0..l2b.tot_len as usize {
            seq[j] = L2b::c2t(l2b.get0(i as u64));
            j += 1;
        }
        for i in 0..l2b.tot_len as usize {
            seq[j] = L2b::g2a(l2b.get0(i as u64));
            j += 1;
        }
        if both_strand {
            for i in (0..l2b.tot_len as usize).rev() {
                seq[j] = 3 - L2b::g2a(l2b.get0(i as u64));
                j += 1;
            }
            for i in (0..l2b.tot_len as usize).rev() {
                seq[j] = 3 - L2b::c2t(l2b.get0(i as u64));
                j += 1;
            }
        }
    } else {
        for i in 0..l2b.tot_len as usize {
            seq[j] = l2b.get0(i as u64);
            j += 1;
        }
        if both_strand {
            for i in (0..l2b.tot_len as usize).rev() {
                seq[j] = 3 - l2b.get0(i as u64);
                j += 1;
            }
        }
    }
    debug_assert_eq!(j, len as usize);

    let n_ssa = ((len + (1u64 << sa_bit)) >> sa_bit) as usize;
    let mut ssa = vec![0u64; n_ssa];
    let threads_i = threads.max(1) as i64;

    let primary = if use_int32 {
        let mut a = vec![0i32; len as usize + FS + 1];
        let ret = if threads <= 1 {
            libsais_rs::libsais(&seq, &mut a[1..], FS as i32, None)
        } else {
            libsais_rs::libsais_omp(&seq, &mut a[1..], FS as i32, None, threads_i as i32)
        };
        if ret != 0 {
            bail!("libsais 构建后缀数组失败（ret={ret}）");
        }
        a[0] = len as i32; // libsais 不写 a[0]，其值恒为 len
        sa_to_bwt_i32(&mut a, &mut seq, len, sa_bit, &mut ssa)
    } else {
        let mut a = vec![0i64; len as usize + FS + 1];
        let ret = if threads <= 1 {
            libsais_rs::libsais64::libsais64(&seq, &mut a[1..], FS as i64, None)
        } else {
            libsais_rs::libsais64::libsais64_omp(&seq, &mut a[1..], FS as i64, None, threads_i)
        };
        if ret != 0 {
            bail!("libsais64 构建后缀数组失败（ret={ret}）");
        }
        a[0] = len as i64;
        sa_to_bwt_i64(&mut a, &mut seq, len, sa_bit, &mut ssa)
    };

    let mut bwt = Bwt::from_raw(&seq, len, primary);
    bwt.sa_bit = sa_bit;
    bwt.sa = ssa;
    Ok(bwt)
}

/// SA → BWT：采样 SA、定位 primary、原地把 SA 转成 BWT 串并移除 $（对应 `sa_to_bwt`）。
fn sa_to_bwt_i32(a: &mut [i32], seq: &mut [u8], len: u64, sa_bit: u32, ssa: &mut [u64]) -> u64 {
    let mask = (1u64 << sa_bit) - 1;
    let mut primary = u64::MAX;
    for i in 0..=len as usize {
        let v = a[i];
        if i as u64 & mask == 0 {
            ssa[i >> sa_bit] = v as u64;
        }
        if v == 0 {
            primary = i as u64;
        } else {
            a[i] = seq[v as usize - 1] as i32;
        }
    }
    ssa[0] = u64::MAX;
    let p = primary as usize;
    for i in 0..p {
        seq[i] = a[i] as u8;
    }
    for i in p..len as usize {
        seq[i] = a[i + 1] as u8;
    }
    primary
}

fn sa_to_bwt_i64(a: &mut [i64], seq: &mut [u8], len: u64, sa_bit: u32, ssa: &mut [u64]) -> u64 {
    let mask = (1u64 << sa_bit) - 1;
    let mut primary = u64::MAX;
    for i in 0..=len as usize {
        let v = a[i];
        if i as u64 & mask == 0 {
            ssa[i >> sa_bit] = v as u64;
        }
        if v == 0 {
            primary = i as u64;
        } else {
            a[i] = seq[v as usize - 1] as i64;
        }
    }
    ssa[0] = u64::MAX;
    let p = primary as usize;
    for i in 0..p {
        seq[i] = a[i] as u8;
    }
    for i in p..len as usize {
        seq[i] = a[i + 1] as u8;
    }
    primary
}

/// 索引加载（对应 `mb_idx_load` 的非 meth 路径）。
pub struct Idx {
    pub l2b: L2b,
    pub bwt: Bwt,
}

impl Idx {
    pub fn load(prefix: &str) -> Result<Idx> {
        let l2b = L2b::load(format!("{prefix}.l2b"))?;
        let mut bwt = Bwt::load(format!("{prefix}.mbw"))?;
        bwt.cache(10); // 对照 mb_idx_load：硬编码 10-mer
        Ok(Idx { l2b, bwt })
    }
}
