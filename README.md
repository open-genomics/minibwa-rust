# minibwa-rust

[minibwa](https://github.com/lh3/minibwa) 的 Rust 重写。minibwa 是 bwa-mem 的继任者：
用 BWT 索引基因组、查找变长种子（SMEM），再以 minimap2 风格的链化与 DP 比对完成短读映射，
比 bwa-mem 快 3 倍以上、比 bwa-mem2 快 2 倍，精度相当。

本项目目标：

- 索引文件与 C 版 minibwa **字节级兼容**（`.l2b` 参考序列 + `.mbw` BWT/采样 SA）
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
| M6 | 配对端：insert size 统计、配对、rescue | 未开始 |
| M7 | 自适应参数全量、对比脚本、CI、文档 | 未开始 |

M5 完成后，本项目已支持单端比对（`map`），并在 chrM 测试数据上与 C 版
输出逐字段一致（含 SA:Z 与全部可选标签）。配对端（M6）尚未完成；
生产环境请使用 [C 版 minibwa](https://github.com/lh3/minibwa)。

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

## 许可

MIT
