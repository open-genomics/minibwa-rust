# 变更记录：PE 双文件读取传播解析错误而非当作记录数不一致

- 日期：2026-11-08
- 变更类型：CLI 错误处理修复

## 内容
- `main.rs` 的 `read_batch_pe` 原先用单个 `match (iter1.next(), iter2.next())`
  只处理 `(Some(Ok), Some(Ok))` 与 `(None, None)`，其余组合（包括任一端
  `Some(Err(..))`）都落入 `_` 分支，打印
  `query files have different number of records; extra records skipped.`
  后 `break`。
- 结果：R1/R2 的解析错误（如 FASTQ qual 长度与 seq 不一致、gzip I/O 错误）
  被误报为记录数不同并丢弃，`run_map` 收到空 batch 后正常结束，进程退出码 0，
  输出只剩 SAM header。
- 修复：任一端返回 `Err` 时立即 `return Err(e.into())`；只有“一端 EOF、
  另一端仍有正常记录”才打印记录数不一致警告并 `break`。

## 影响
- 非法 FASTQ 输入不再静默成功：CLI 输出具体解析错误（如
  `Error: qual 长度与 seq 不一致`）并以非 0 退出码结束。
- 合法输入的读取/配对路径逐字节不变；单端路径 `read_batch` 原本就传播错误。

## 验证
- 新增 `src/main.rs::tests::read_batch_pe_propagates_read_errors`：R1 为
  `@a\nACGT\n+\nIII\n`（qual 比 seq 短），R2 合法；修复前返回 Ok，
  修复后返回 Err。
- CLI 复核（/tmp 夹具 + 现有 /tmp/mini 索引）：
  `minibwa-rust map -t 1 /tmp/mini /tmp/pe-bad-r1.fq /tmp/pe-good-r2.fq`
  修复后 exit=1、stderr 为 `Error: qual 长度与 seq 不一致`；合法的
  `@a/1` + `@a/2` 配对仍 exit=0 且正常输出两条 SAM 记录。
- `cargo test -j 4`：88 passed / 0 failed；`cargo clippy --all-targets --
  -D warnings`、`cargo fmt --all -- --check` 通过。
- 本地缺 `../minibwa`，16 个 c_parity 测试全部 `[PARITY SKIP]`；本修复只影响
  错误路径，合法输入的 C 对照行为不变。

## 路径
- 代码：`src/main.rs`
- 测试：`src/main.rs`（tests 模块）
- 本记录：`changelog/2026-11-08-pe-read-error-propagation.md`
