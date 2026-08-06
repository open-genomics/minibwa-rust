//! l2bit 参考基因组格式（`.l2b`），对应 C 版 `l2bit.c`。

use anyhow::Result;

/// 从 `.l2b` 还原参考基因组并打印到 stdout（对应 C 版 `getref`）。
pub fn getref(_path: &str) -> Result<()> {
    anyhow::bail!("l2bit 模块尚未实现（M2 里程碑）")
}
