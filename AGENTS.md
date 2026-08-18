# AGENTS.md

## Project Truth

minibwa-rust 是 C 版 [minibwa](https://github.com/lh3/minibwa)（bwa-mem 继任者）的 Rust 重写。
参考 C 源码位于 `/home/shane/github/open-genomics/minibwa`（跟踪 pin 住的 commit，勿改）。

核心承诺：**索引文件（.l2b/.mbw）与 C 版字节级兼容**，算法行为逐模块对照。

## Source Of Truth

- 参数默认值：`src/opt.rs`（对照 C 版 options.c 的 mb_opt_reset/mb_opt_preset/mb_opt_adap）
- 索引格式：`src/l2bit.rs`（l2bit.c）、`src/bwt.rs`（bwt.c）
- CLI：`src/main.rs`
- C 参考 commit：仓库根目录 `C_MINIBWA_COMMIT`

## Conventions

- 纯 Rust，`unsafe_code = "forbid"`；后缀数组用 libsais-rs
- 字段命名与 C 版保持一致（sa_bit、min_len、max_sub_occ…）以便逐行对照
- 提交：语义化原子提交，中文说明；每个提交独立可构建、cargo test 全绿
- 验证：凡与 C 版有交互的格式/算法，必须有与 C 二进制的输出对照测试
- 文档只保留中文，避免过度工程化

## Verification

```bash
# 准备 C 参考（commit 见 C_MINIBWA_COMMIT）
git -C ../minibwa checkout "$(grep -v '^#' C_MINIBWA_COMMIT | tr -d '[:space:]')"
make -C ../minibwa

cargo test
cargo clippy --all-targets -- -D warnings
./scripts/compare.sh <ref.fa> <reads.fq>   # C minibwa map vs Rust map
```

本地缺少 `../minibwa/minibwa` 时，`c_parity` 打印 `[PARITY SKIP]` 并跳过。
CI 设置 `MINIBWA_REQUIRE_PARITY=1`，会编译 pin 住的 C 二进制；缺依赖即失败。
