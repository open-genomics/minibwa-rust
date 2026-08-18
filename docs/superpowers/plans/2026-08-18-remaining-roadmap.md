# 后续路线图（P2 剩余 + P3）

> **For agentic workers:** P2a（10-mer cache）已完成。本文件记录剩余项与推荐切片。SIMD / 流水线 / `--meth` 均需用户确认后再写代码。

**Goal:** 在 cache 对齐之后，把剩余缺口分成可独立交付的切片，避免把性能优化和 BS-seq 搅在一起。

**Architecture:** 短读路径算法与 C 对照已闭环。剩下三类工作互不依赖：SIMD 政策、I/O 流水线、`--meth` 子系统。优先功能缺口（`--meth`），性能项单独决策。

**Tech Stack:** 现有纯 Rust crate；`unsafe_code = "forbid"`；C 参考 pin 在 `C_MINIBWA_COMMIT`。

**Spec:** 本对话的完成度分析；C `index.c` / `l2bit.c` / `map-algo.c` / `map-main.c` / `ksw2_extd2_sse.c`。

## 已完成（本会话及之前）

| 阶段 | 内容 | 证据 |
|------|------|------|
| P0 | pin C、CI parity、`-x sr`、SE SAM/PAF、`compare.sh` | `tests/c_parity.rs` |
| P1 | `--eqx` / `-b` / `--hic` / `--adap` / `--pe` / 打分 / RG | 5 个 map 对照测试 |
| P2a | `back()` 游标 + `Idx::load` 启用 `cache(10)` | `smem_with_cache_matches_without`；`MINIBWA_REQUIRE_PARITY=1 cargo test --all-targets --all-features` 64+15 绿 |

## 本计划不做（除非用户明确点名）

- GPL `index -l` / bwtgen
- 改 C 参考仓库

---

### 切片 A — SIMD ksw2（阻塞）

C 用 SSE/NEON 条带 DP（`ksw2_extd2_sse.c`）。Rust `src/ksw2.rs` 是标量、分数与 CIGAR 已与 chrM 对照一致。

**阻塞：** `[lints.rust] unsafe_code = "forbid"`。稳定版 Rust 没有安全的 `std::simd`；依赖 crate 内部 unsafe 不违反本 crate 的 forbid，但条带 DP 要对齐 C 的数值行为，工作量大、容易破坏 parity。

**建议：** 暂不移植。若要做，先改政策（允许 `unsafe` 于 `ksw2` 模块，或引入安全 SIMD crate），并加独立 DP 对照 harness，不要先改 map 路径。

### 切片 B — 读/比对/写流水线（性能，非正确性）

C `map-main.c`：`kt_pipeline` 三步（读 → `kt_for` 比对 → 写），读写下一批与比对重叠。

Rust `run_map`：按 mini-batch 串行读/写，batch 内 `rayon` 并行 `map_one`。chrM 对照已与 C SAM 一致。缺的是 **I/O 与计算重叠**，不是并行度本身。

**建议：** 没有大样本吞吐基准前不改。改动集中在 `src/main.rs`，必须保持 SAM 记录顺序，否则 parity 会碎。

### 切片 C — `--meth` BS-seq（下一功能缺口，架构级）

C 行为要点：

1. `index --meth`：照常写 `.l2b` + `.mbw`，再对 C2T∥G2A（各正反向）建 `.meth.mbw`（`mb_bwt_libsais(..., is_meth=1)`）。
2. `map --meth`：加载 `.meth.mbw`；R1 查询 C→T、R2 G→A；锚点 `l2b_intv2cid_meth`；打分矩阵把 C-T / G-A 当匹配（`gen_nt4_mat` 的 `mt` 已在 Rust 里）。
3. SAM：XR/XG/XM 等甲基化标签。

Rust 现状：`index::build` / `Idx::load` / `mapcore` / `seed` 全部标了「非 meth」。`ksw2::gen_nt4_mat` 已有 `mt` 参数，未接入映射。

**推荐交付顺序（每步可独立对照 C）：**

1. `index --meth`：`.meth.mbw` 与 C 字节级一致。**已完成**（`index_parity_meth`）。
2. `Idx::load(..., is_meth)` 读 `.meth.mbw` + `cache(10)`。
3. 查询转换 + 两轮 SMEM/锚点 meth 路径，chrM 或 toy FASTA 上 `map --meth` SAM 对照。
4. PE rescue 的 meth 分支 + XR/XG/XM。

## 推荐下一步

用户确认后从 **切片 C.1：`index --meth` 字节对照** 开始。A/B 保持文档状态，不写代码。
