//! 索引构建：FASTA → `.l2b` + `.mbw`，对应 C 版 `index.c`。

use anyhow::Result;

/// 为参考基因组 FASTA 建索引（对应 C 版 `main_index` 的 libsais 路径）。
pub fn build(_fasta: &str, _prefix: &str, _seed: u64, _sa_bit: u32, _threads: usize) -> Result<()> {
    anyhow::bail!("index 模块尚未实现（M2 里程碑）")
}
