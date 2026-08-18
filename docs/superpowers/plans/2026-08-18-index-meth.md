# C.1：`index --meth` 字节对照 Implementation Plan

> **For agentic workers:** 本会话按用户确认：只做切片 C.1。写计划后立即 inline 执行。不做 map `--meth`、`Idx::load` meth、SIMD、流水线。

**Goal:** `minibwa-rust index --meth` 在写出常规 `.l2b`/`.mbw` 的同时，再建与 C 版字节一致的 `.meth.mbw`。

**Architecture:** C 非 `-l` 路径：`mb_bwt_libsais(l2b, sa_bit, both_strand=1, is_meth=0)` 写 `.mbw`，再 `is_meth=1` 写 `.meth.mbw`。meth 文本是 C2T 正向 ∥ G2A 正向 ∥ G2A 反向互补 ∥ C2T 反向互补。Rust 只给 `bwt_libsais` 加 `is_meth`，序列填充对照 `index.c`。

**Tech Stack:** 现有 `src/index.rs` / `src/l2bit.rs` / `src/main.rs`，`tests/c_parity.rs`。

**Spec:** 对话中的切片 C.1；C `index.c` `mb_bwt_libsais` / `main_index`。

## Global Constraints

- `unsafe_code = "forbid"`
- 字段命名与 C 一致
- 有 C 交互的格式必须有对照测试
- 文档中文
- 不提交除非用户要求

## 本计划不做

- `Idx::load` 读 `.meth.mbw`
- `map --meth` / 查询转换 / 锚点 / SAM 标签
- GPL `index -l`

---

### Task 1: 对照测试先红

**Files:**
- Modify: `tests/c_parity.rs`

- [x] 新增 `index_parity_meth`：C/Rust 均 `index --meth -t1`，比较 `.l2b`、`.mbw`、`.meth.mbw`
- [x] 跑测试，确认因缺少 `--meth` 或未写 `.meth.mbw` 而失败

### Task 2: 实现 meth 文本 + CLI

**Files:**
- Modify: `src/l2bit.rs`（`c2t`/`g2a`，对照 `l2b_c2t`/`l2b_g2a`）
- Modify: `src/index.rs`（`bwt_libsais(..., is_meth)`；`build` 在 meth 时再写 `.meth.mbw`）
- Modify: `src/main.rs`（`index --meth`）

- [x] 最小改动让 `index_parity_meth` 转绿
- [x] 既有 `index_parity_*` 仍绿（无 `--meth` 不写 `.meth.mbw`）

### Task 3: 验证与文档

- [x] `MINIBWA_REQUIRE_PARITY=1 cargo test --all-targets --all-features`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] README：`index --meth` 生成 `.meth.mbw`（映射仍未做）
