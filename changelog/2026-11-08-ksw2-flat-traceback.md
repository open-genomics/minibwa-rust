# 2026-11-08 ksw2 回溯矩阵改为扁平缓冲

## 问题

`ksw2::extz` / `extd2` 的回溯矩阵按反斜角逐行保存：每个反斜角把
`drow: Vec<u8>` push 进 `p: Vec<Vec<u8>>`。一次比对的反斜角数约为
`(qlen + tlen) / 2`，长读（数 kb）动辄上千次小块分配，回收时同样上千次释放；
矩阵内容本身也不连续，回溯时每步都要多一次间接寻址。

## 改动

- `p` 改为单个 `Vec<u8>`，每行 `extend_from_slice` 追加到末尾；新增
  `p_off: Vec<usize>` 记录每行起始偏移（行长度仍由 `off`/`off_end` 决定，
  与原来每行 `en - st + 1` 一致）。
- `backtrack` 签名改为 `p: &[u8], p_off: &[usize], ...`，取字节改为
  `p[p_off[r] + (i - off[r]) as usize]`。
- 三个调用点（正常回溯、extz-only 回溯、best-hit 回溯）同步传 `p_off`。

分配次数从「每个反斜角一次」降到「容量翻倍时一次」，内存布局连续；
`drow` 仍按行构造（`vec![0u8; 行宽]`），未改动 DP 数值路径，因此 CIGAR
输出不受影响。

## 验证

- `cargo test --quiet`：68 lib + 1 bin + 16 C 对拍 + 9 损坏索引 + 1 集成全部通过。
  C 对拍（`tests/c_parity`）逐条比较与 ksw2 C 实现的 CIGAR/得分，覆盖本次改动。
- `cargo clippy --all-targets -- -D warnings`、`cargo fmt --all -- --check` 通过。

## 未做

- `drow` 本身仍每行分配一次；进一步可复用一块临时行缓冲（`clear()` 后
  `resize`），但收益只剩一半且会引入「行宽变化时的容量管理」，暂缓。
