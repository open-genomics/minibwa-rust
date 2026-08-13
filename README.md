# minibwa-rust

[minibwa](https://github.com/lh3/minibwa) 的 Rust 重写。minibwa 是 bwa-mem 的继任者：
用 BWT 索引基因组、查找变长种子（SMEM），再以 minimap2 风格的链化与 DP 比对完成短读映射，
比 bwa-mem 快 3 倍以上、比 bwa-mem2 快 2 倍，精度相当。

本项目目标：

- 索引文件与 C 版 minibwa **字节级兼容**（`.l2b` 参考序列 + `.mbw` BWT/采样 SA）

> **Parity 测试前提**：`tests/c_parity.rs` 需要已编译的 C 版 minibwa 二进制
> （`../minibwa/minibwa`）。该二进制不存在时，parity 测试会跳过并打印提示，
> **不视为 C parity 已验证**。当前 CI 未构建 C 参考二进制，因此 CI 通过
> 仅表示 Rust 自身测试通过，不证明与 C 版的字节级一致性。
- 算法逐模块对照移植，以 C 版输出为基准做一致性验证
- 纯 Rust 实现（后缀数组构建使用 libsais 的纯 Rust 移植 `libsais-rs`），无 C 依赖

## 当前状态

| 里程碑 | 内容 | 状态 |
|--------|------|------|
| M1 | 仓库骨架：CLI、参数体系（`mb_opt_t`）、DNA 工具 | ✅ |
| M2 | l2bit + BWT + 索引构建（`index` / `getref`） | ✅ |
| M3 | 两轮 SMEM 种子 + kmer 缓存 | ✅ |
| M4 | 链化 + hit 处理 + mapQ | ✅ |
| M5 | ksw2 风格比对 + CIGAR + SAM/PAF（单端一致） | ✅ |
| M6 | 配对端：insert size 统计、配对、rescue | ✅ |
| M7 | 自适应参数全量、对比脚本、CI、文档 | 未开始 |

M6 完成后，本项目已支持配对端比对（`map` 与 `mem` 双命令），并在 chrM 测试
数据上与 C 版输出逐字节一致：单端、配对端（含 pestat 数据估计、proper pair
配对、mate rescue、PE mapQ）、`-P` 跳过配对、单文件交错输入、`mem` 的
MD:Z / n2 / XA 标签。

## 构建与使用

```bash
cargo build --release

# 建索引（与 C 版 minibwa index 兼容）
minibwa-rust index ref.fa prefix        # 生成 prefix.l2b + prefix.mbw

# 从 .l2b 还原参考序列（调试）
minibwa-rust getref prefix.l2b

# 比对
minibwa-rust map prefix read1.fq read2.fq > aln.sam
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
| `src/dna.rs` | kommon.c（部分） | 编码表、splitmix64 |

## 测试

```bash
cargo test
```

- M2：索引字节级兼容性验证（`tests/` 集成测试）
- M5：chrM 单端 SAM/PAF 与 C 版 minibwa 逐字段一致（含 SA:Z、NM/AS/ms/md、cg:Z）
- M6：chrM 配对端 SAM 与 C 版逐行一致（`tests/c_parity.rs` 的 `pe_sam_parity` /
  `mem_pe_parity`，覆盖 pestat、配对、rescue、`-P`、`-5`、`-I`、`-S`）

## 许可

MIT
