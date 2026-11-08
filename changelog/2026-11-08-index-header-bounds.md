# 变更记录：损坏索引头不再驱动巨量分配

- 日期：2026-11-08
- 变更类型：健壮性修复（纯 Rust，索引格式不变）

## 内容
- `Bwt::load` 在分配 `data` 前用文件剩余字节校验 `data_len(seq_len)`：
  `seq_len` 来自 `.mbw` 头部（`l2[4]`），旧实现直接
  `vec![0u64; data_len(seq_len)]`，把 `l2[4]` 改成 `1<<50` 的几 KB 文件会
  触发数百 TB 分配并 abort。
- 读取 `n_sa` 后同样用剩余字节校验；并在计算
  `(seq_len + (1 << sa_bit)) >> sa_bit` 前拒绝 `sa_bit >= 64`（移位在 Rust
  中会 panic/回绕）。
- `L2b::from_bytes` 把头部 7 个计数全部按固定记录尺寸求和
  （`n_ctg×8 + (n_ambi+n_mask)×16 + n_pac×8 + len_name + len_comm`），
  与文件剩余字节核对后才 `Vec::with_capacity` / `vec![..; n]`。
  另加：
  - `n_pac == ceil(tot_len/32)` 一致性检查（防止后续 `get0` 越界）；
  - contig 长度累加改用 `checked_add`。
- 新增 `tests/test_corrupt_index.rs` 5 条回归：`.mbw` 巨大 `seq_len`、
  非法 `sa_bit`；`.l2b` 巨大 `n_ctg`、巨大 `n_pac`、`n_pac` 不一致。

## 影响
- 损坏/伪造的 `.mbw` / `.l2b` 在读取数据前返回可诊断的 `anyhow` 错误，
  而不是 `memory allocation failed` abort 或 `capacity overflow` panic。
- 合法索引的解析路径与文件格式不变；`n_pac` 检查按 C 版写入约定
  `(tot_len + 31) / 32` 执行。
- 大文件本身仍需要与其内容相称的内存；本次只消除“小文件声明大尺寸”的
  放大效应。

## 验证
- `cargo test`：86 passed / 0 failed（含 `c_parity` 在缺少 C 二进制时按
  约定跳过）；新增 5 条测试全部通过。
- `cargo fmt --all -- --check`、`cargo clippy --all-targets -- -D warnings`：
  通过。
- 未用 C 版 minibwa 生成真实索引做交叉加载（本机无 `../minibwa` 构建）；
  格式语义未改动。

## 路径
- 代码：`src/bwt.rs`、`src/l2bit.rs`
- 测试：`tests/test_corrupt_index.rs`
- 本记录：`changelog/2026-11-08-index-header-bounds.md`
