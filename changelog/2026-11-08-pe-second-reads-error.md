# 变更记录：PE 第二个 reads 文件打开失败返回错误而非 panic

- 日期：2026-11-08
- 变更类型：CLI 错误处理修复

## 内容
- `main.rs` 的 PE 路径原本写：

  ```rust
  let mut it2 = reads2.map(|p| SeqIter::new(open_fasta(p).unwrap()));
  ```

  第二个 reads 文件不存在/不可读时 `unwrap` 直接 panic（release 下
  `panic=abort` 为 SIGABRT），与同函数其他 I/O 路径的 `anyhow` 错误不一致。
- 改为 `match` + `?` 传播：

  ```rust
  let mut it2 = match reads2 {
      Some(p) => Some(SeqIter::new(open_fasta(p)?)),
      None => None,
  };
  ```

## 影响
- 用户传入不存在的 R2 文件时，CLI 输出正常错误并以一致的退出码结束。
- 单端路径（`reads2 == None`）与成功路径行为不变。

## 验证
- `cargo test`：86 passed / 0 failed；`cargo clippy --all-targets -- -D warnings`、
  `cargo fmt --check` 通过。
- 未单独启动 CLI 复现“R2 不存在”场景（需要完整索引夹具）；改动只把
  `unwrap` 换成 `?`，编译期已覆盖类型/错误传播。

## 路径
- 代码：`src/main.rs`
- 本记录：`changelog/2026-11-08-pe-second-reads-error.md`
