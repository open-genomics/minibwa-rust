# 2026-11-08 `.mbw` 头部语义校验

`Bwt::load` 在原有长度/剩余字节检查之外，补齐文件头语义校验：

- `seq_len == 0` 拒绝；
- `primary >= seq_len` 拒绝；
- `L2` 必须非递减；
- `sa_bit` 有效时 `n_sa` 必须等于推导值（含 `n_sa == 0` 的拒绝，避免 map 阶段
  `cache(10)` 越界）；`sa_bit == u32::MAX` 时 `n_sa` 必须为 0。

`tests/test_corrupt_index.rs` 新增 4 条回归：零 SA、primary 越界、L2 非单调、
非法 sa_bit（原有用例按新前置条件调整）。

## 验证

`cargo test -j 4`：65 lib + 1 bin + 16 c_parity + 8 corrupt_index + doctest 全部通过；
`cargo clippy --all-targets -- -D warnings` 通过。`c_parity` 因本地无 `../minibwa`
仍为 SKIP。
