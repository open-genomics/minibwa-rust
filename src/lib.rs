//! # minibwa-rust
//!
//! [minibwa](https://github.com/lh3/minibwa) 的 Rust 重写。
//!
//! minibwa 是 bwa-mem 的继任者：用 BWT 索引基因组、像 bwa-mem 一样查找
//! 变长种子（SMEM），再使用 minimap2 的链化与 SIMD 比对算法。
//! 本项目追求与 C 版 minibwa 的索引文件（`.l2b` / `.mbw`）字节级兼容。
//!
//! ## 模块说明
//!
//! - [`dna`] — DNA 编码表、反向互补、splitmix64 随机数
//! - [`opt`] — 比对参数（对应 C 版 `mb_opt_t`），含 sr/lr/adap 预设
//! - [`io`] — FASTA / FASTQ 解析
//! - [`l2bit`] — l2bit 参考基因组格式（`.l2b`）
//! - [`bwt`] — BWT + 采样后缀数组（`.mbw`），Rank / SMEM / SA 查询
//! - [`index`] — 索引构建（FASTA → `.l2b` + `.mbw`）

pub mod bwt;
pub mod dna;
pub mod index;
pub mod io;
pub mod l2bit;
pub mod opt;
pub mod chain;
pub mod format;
pub mod mapcore;
pub mod hits;
pub mod seed;
