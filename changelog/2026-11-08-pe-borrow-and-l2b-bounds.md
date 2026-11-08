# 2026-11-08 PE 配对借用序列、l2b 越界防御与依赖升级

## PE 配对不再深拷贝 read

`main.rs` 的逐 fragment 配对循环原先为每个 fragment 克隆两端序列
（`batch.seqs[..].seq.clone()`），只为满足 `pair()` 的
`&[Vec<u8>; 2]` 参数。`pair()` / `matesw()` 只读序列，参数改为
`&[&[u8]; 2]`，调用方直接借用 `batch.seqs`，每个 fragment 省掉两次
`String` 深拷贝。

## `.l2b` 名字/注释串越界

头部 `len_name` / `len_comm` 小于 contig 数时，旧实现直接
`names[p_name..]` 与 `comms[p_comm]` 索引，恶意索引可让 `map/getref`
加载即 panic。现在用 `get(..)` / `get(p_comm)` 检查并返回可诊断格式错误；
`tests/test_corrupt_index.rs` 新增 `len_comm=0, n_ctg=1` 回归。

## 依赖 advisory

`cargo update -p crossbeam-epoch -p anyhow`：

- `crossbeam-epoch` 0.9.18 → 0.9.21（RUSTSEC-2026-0204）
- `anyhow` 1.0.102 → 1.0.104（RUSTSEC-2026-0190）

## 验证

- `cargo test --all-targets`：68 lib + 1 bin + 16 c_parity + 9 corrupt_index
  全部通过；
- `cargo clippy --all-targets -- -D warnings`、`cargo fmt --all -- --check` 通过。
