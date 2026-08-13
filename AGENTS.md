# AGENTS.md

## Project Truth

minibwa-rust 是 C 版 [minibwa](https://github.com/lh3/minibwa)（bwa-mem 继任者）的 Rust 重写。
参考 C 源码位于 `/home/shane/github/open-genomics/minibwa`（跟踪 upstream，勿改）。

核心承诺：**索引文件（.l2b/.mbw）与 C 版字节级兼容**，算法行为逐模块对照。

## Source Of Truth

- 参数默认值：`src/opt.rs`（对照 C 版 options.c 的 mb_opt_reset/mb_opt_preset/mb_opt_adap）
- 索引格式：`src/l2bit.rs`（l2bit.c）、`src/bwt.rs`（bwt.c）
- CLI：`src/main.rs`

## Conventions

- 纯 Rust，`unsafe_code = "forbid"`；后缀数组用 libsais-rs
- 字段命名与 C 版保持一致（sa_bit、min_len、max_sub_occ…）以便逐行对照
- 提交：语义化原子提交，中文说明；每个提交独立可构建、cargo test 全绿
- 验证：凡与 C 版有交互的格式/算法，必须有与 C 二进制的输出对照测试
- 文档只保留中文，避免过度工程化

## Verification

```bash
cargo test
cargo clippy --all-targets -- -D warnings
# C 版对照（需要 ../minibwa/minibwa 二进制；不存在时 parity 测试跳过）
# 当前 CI 未构建 C 参考二进制，CI 通过 ≠ C parity 已验证
./scripts/compare.sh <ref.fa> <reads.fq>   # 注意：脚本使用 bwa（非 minibwa），
                                           # 且 Rust CLI 接口已漂移，待修复
