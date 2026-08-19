# minibwa-rust

[minibwa](https://github.com/lh3/minibwa) 的 Rust 重写。minibwa 是 bwa-mem 的继任者：
用 BWT 索引基因组、查找变长种子（SMEM），再以 minimap2 风格的链化与 DP 比对完成短读映射，
比 bwa-mem 快 3 倍以上、比 bwa-mem2 快 2 倍，精度相当。

本项目目标：

- 索引文件与 C 版 minibwa **字节级兼容**（`.l2b` 参考序列 + `.mbw` BWT/采样 SA；`--meth` 时另有 `.meth.mbw`）
- 算法逐模块对照移植，以 **已 pin 的 C 版 commit** 输出为基准做一致性验证
- 纯 Rust 实现（后缀数组构建使用 libsais 的纯 Rust 移植 `libsais-rs`），无 C 依赖
- **本项目已功能完备（标准 DNA 比对），性能优化与 BS-seq 映射不在范围内**

> **C 参考**：`https://github.com/lh3/minibwa`，commit 写在仓库根目录 `C_MINIBWA_COMMIT`。
> 本地路径 `../minibwa`（相对本仓库根）。构建：`make -C ../minibwa`。
>
> **Parity**：`tests/c_parity.rs` 需要 `../minibwa/minibwa` 与 `../minibwa/test/chrM-*.gz`。
> 本地未准备时打印 `[PARITY SKIP]` 并跳过（**不视为已验证**）。
> CI 设置 `MINIBWA_REQUIRE_PARITY=1` 并编译 pin 住的 C 二进制，缺依赖即失败。

## 当前状态

| 里程碑 | 内容 | 状态 |
|--------|------|------|
| M1 | 仓库骨架：CLI、参数体系（`mb_opt_t`）、DNA 工具 | ✅ |
| M2 | l2bit + BWT + 索引构建（`index` / `getref`） | ✅ |
| M3 | 两轮 SMEM 种子 + kmer 缓存 | ✅ |
| M4 | 链化 + hit 处理 + mapQ | ✅ |
| M5 | ksw2 风格比对 + CIGAR + SAM/PAF（单端一致） | ✅ |
| M6 | 配对端：insert size 统计、配对、rescue | ✅ |
| M7 | 验证闭环（pin C、CI parity、对照脚本、文档） | ✅ |
| M8 | 常用 map CLI（`--eqx`/`-b`/`--hic`/`--adap`/`--pe`） | ✅ |
| M9 | `Idx::load` 启用 10-mer SMEM cache（与 C `mb_bwt_cache(bwt, 10)` 对齐） | ✅ |
| M10 | `index --meth` 写出与 C 字节一致的 `.meth.mbw` | ✅ |
| — | 最终收尾：compare.sh 修复、文档定型 | ✅ |

## 已知限制

以下功能在 C 版 minibwa 中存在，但本 Rust 重写因政策或优先级原因不实现：

| 功能 | 原因 |
|------|------|
| SIMD ksw2（SSE/NEON 条带 DP） | 与 `unsafe_code = "forbid"` 冲突；标量版已与 C 对照一致 |
| 读/比对/写流水线 | 纯性能优化，不影响正确性；当前 batch 内 rayon 并行已足够 |
| `map --meth` BS-seq 映射 | `index --meth` 已与 C 字节一致，映射路径未实现 |

M6 完成后，本项目已支持配对端比对（`map` 与 `mem` 双命令），并在 chrM 测试
数据上与 C 版输出对照：单端 SAM/PAF、配对端（含 pestat、proper pair、mate rescue、
PE mapQ）、`-P` 跳过配对、`mem` 的 MD:Z / n2 / XA 标签。

## 构建与使用

```bash
cargo build --release

# 建索引（与 C 版 minibwa index 兼容）
minibwa-rust index ref.fa prefix        # 生成 prefix.l2b + prefix.mbw
minibwa-rust index --meth ref.fa prefix # 额外生成 prefix.meth.mbw（映射尚未接入）

# 从 .l2b 还原参考序列（调试）
minibwa-rust getref prefix.l2b

# 比对
minibwa-rust map prefix read1.fq read2.fq > aln.sam
```

准备 C 对照：

```bash
git clone https://github.com/lh3/minibwa.git ../minibwa
git -C ../minibwa checkout "$(grep -v '^#' C_MINIBWA_COMMIT | tr -d '[:space:]')"
make -C ../minibwa
cargo test --test c_parity
./scripts/compare.sh ref.fa reads.fq
```

## 架构

```text
FASTA ──> l2bit(.l2b) ──> libsais SA ──> BWT + 采样 SA(.mbw)
                                              │
reads ──> 两轮 SMEM 种子 ──> 链化 ──> ksw2 比对 ──> SAM/PAF
```

模块与 C 版对照：

| Rust 模块 | C 版文件 | 职责 |
|-----------|----------|------|
| `src/l2bit.rs` | l2bit.c | 2-bit 参考基因组格式读写 |
| `src/bwt.rs` | bwt.c | BWT Rank/extend/SMEM、采样 SA |
| `src/index.rs` | index.c | 索引构建流水线 |
| `src/opt.rs` | options.c | 参数默认值、sr/lr/adap 预设 |
| `src/cs.rs` | cs.c | cs:Z / ds:Z 标签 |
| `src/dna.rs` | kommon.c（部分） | 编码表、splitmix64 |

## 测试

```bash
cargo test
```

- M2：索引字节级兼容（`index_parity_*` / `getref_parity`）
- M10：`index --meth` 的 `.meth.mbw` 字节对照（`index_parity_meth`）
- M5：chrM 单端 `se_sam_parity` / `se_paf_parity`（去 header 后逐行 / `-f` 逐字节）
- M6：chrM 配对端 `pe_sam_parity` / `mem_pe_parity`（pestat、配对、rescue、`-P`、`-5`、`-I`、`-S`）

## 许可

MIT
