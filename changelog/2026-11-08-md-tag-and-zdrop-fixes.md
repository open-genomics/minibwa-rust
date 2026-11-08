# 2026-11-08 MD:Z 前导 indel、Z-drop 回溯与分裂尾部比对

三处与 C 版 `align.c` 对齐的修复。

## MD:Z 使用 fix_cigar 平移后的序列

C 版 `mb_update_extra` 在调用 `mb_fix_cigar` 后立即执行
`qseq += qshift, tseq += tshift`，后续打分、EQX 与 MD/cs 全部使用平移后的
序列。Rust 版此前额外保留了 `qseq_full`/`tseq_full` 专门给 MD 用，导致
CIGAR 的前导 I/D 被删除后序列仍从原起点计算，MD 标签整体错位
（例如应输出 `3` 却输出 `0^TT3`）。

现在直接使用平移后的 `qseq`/`tseq`，与 C 版一致。新增两个单元测试：
前导 2D、前导 2I 都断言 `MD:Z:3`；另保留一个普通错配用例
（`MD:Z:1C1`）。

## Z-drop 回溯起点

`align1` 在 `ez.zdropped` 后寻找最后一个满足
`a[j].tpos <= ts + ez.max_t` 的锚点时，初值写成了 `i - 2`，跳过了
`i - 1`。C 版从 `i - 1` 开始，现已修正。

## align_skeleton 循环上界

`align_skeleton` 的逐 hit 比对循环把上界固定为初始 `hits.len()`。C 版的
`for (i = 0; i < n_regs; ++i)` 使用会被 `mb_insert_reg` 增加的 `n_regs`，
因此 Z-drop 分裂出的尾部 hit 也会执行 `align1`。Rust 版新插入的尾部 hit
从未做碱基比对，没有 CIGAR/MD，随后会被过滤掉。循环改为
`while i < hits.len()`，插入后下一轮继续处理，与 C 版行为一致。

## 验证

- `cargo test --all-targets`：68 lib + 1 bin + 16 c_parity + 8 corrupt_index 全部通过。
- `cargo clippy --all-targets -- -D warnings`、`cargo fmt --all -- --check` 通过。
- 环境无本地 `../minibwa` 时 `c_parity` 仍按原有逻辑跳过。
