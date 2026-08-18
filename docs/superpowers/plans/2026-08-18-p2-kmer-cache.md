# 10-mer SMEM cache 对齐 Implementation Plan

> **For agentic workers:** 本会话按用户指示：写计划后立即 inline 执行。

**Goal:** 让 `Bwt::cache` + `back()` 与无缓存 SMEM 结果一致，并在 `Idx::load` 启用与 C 版相同的 10-mer 缓存。

**Architecture:** C 的 `mb_bwt_back` 在命中 `pre[]` 后把游标留在 `pos - pre_len`；Rust 用独立变量 `j` 收集 k-mer，命中后 `i` 仍是 `pos - 1`，把已缓存的碱基又 extend 了一遍。只改 `back()` 的游标，并在 `Idx::load` 调用 `cache(10)`（对照 `mb_idx_load`）。

**Tech Stack:** Rust, 现有 `src/bwt.rs` / `src/index.rs`，chrM C parity。

**Spec:** 对话中的 P2 第一项；C `bwt.c` `mb_bwt_back` / `mb_idx_load`。

## Global Constraints

- `unsafe_code = "forbid"`
- 字段命名与 C 一致
- 有 C 交互的算法必须有对照测试
- 文档中文

## 本计划不做

- SIMD ksw2（与 `unsafe_code = "forbid"` 冲突，需单独政策）
- `--meth`
- 读/比对/写流水线
- `mb_seed_intv_batch` 状态机移植（cache 正确后标量 SMEM 已与 C 输出对齐）

---

### Task 1: 回归测试 — 有缓存与无缓存 SMEM 一致

**Files:**
- Modify: `src/bwt.rs` 测试模块

- [x] 写 `smem_with_cache_matches_without`：`min_len` 足够触发 cache（`min_len - 1 >= pre_len`）
- [x] 跑测试，确认因游标错误而失败
- [x] 修 `back()`：cache 命中后 `i = pos - pre_len`（与 C 的 `for (i = pos; l < pre_len; --i)` 结束值一致）
- [x] 测试转绿

### Task 2: Idx::load 启用 cache(10)

**Files:**
- Modify: `src/index.rs`

- [x] `Bwt::load` 后 `bwt.cache(10)`，对照 `mb_idx_load`
- [x] `MINIBWA_REQUIRE_PARITY=1 cargo test --test c_parity` 全绿

### Task 3: 文档

**Files:**
- Modify: `README.md`

- [x] 里程碑表：10-mer cache 已在 load 路径启用；SIMD / meth 仍未做
