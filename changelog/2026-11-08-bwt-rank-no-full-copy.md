# 变更记录：BWT rank 不再复制整段索引

- 日期：2026-11-08
- 变更类型：性能修复（纯 Rust，无格式变化）

## 内容
- `Bwt::rank1a` 不再对 `&self.data[block + 4..]` 做
  `flat_map(...).collect::<Vec<u32>>()`。原实现每次 rank 调用都会把当前
  128 bp 块之后的**整个 BWT 数据段**展开并复制成 `Vec<u32>`，复杂度
  O(N)，N 为索引长度；实际只用到块内 8 个 u32。
- 新增 `#[inline] fn word32(base, i)`，按 `i >> 1` / `i & 1` 直接读取
  块内第 i 个 u32（小端），`rank1a` 与 `rank2a` 都改为按需读取。
- `rank2a` 原先为 4 个 u64 建一个临时 `Vec<u32>`，同样改为直接下标读取。

## 影响
- rank 的每次调用从 O(N) 分配+复制降为 O(1) 读取，SMEM/seed 热路径不再
  随参考基因组长度线性变慢。
- 索引文件（`.l2b` / `.mbw`）与 C 版字节级兼容性、rank 数值语义均不变。
- `bwt.rs` 的公开 API 不变，新增的 `word32` 为私有辅助函数。

## 验证
- `cargo test`：81 passed / 0 failed（含 c_parity 在缺少 C 二进制时按约定跳过）。
- `cargo fmt --check`、`cargo clippy --all-targets -- -D warnings`：通过。
- 同一 1 Mbp 参考与 20k reads（3 条模板循环），`map -t 1`、release、
  各跑 3 次：

  | 版本 | wall time | 最大 RSS |
  |------|-----------|----------|
  | HEAD（旧实现） | 9.51 / 9.60 / 9.66 s | ~68 MB |
  | 本次（word32） | 0.09 / 0.09 / 0.09 s | ~67 MB |

  旧实现约 9.6 s，新实现约 0.09 s，约 106×；200k reads 新实现 0.70 s。
  两个版本的 SAM 输出逐字节相同（20003 行）。计时含进程启动与索引加载，
  不含索引构建；未在真实基因组规模上重复该基准。
- 未运行 `scripts/compare.sh` 的 C 版对照；本次不改变任何输出语义，
  SAM 新旧一致性已覆盖。

## 路径
- 代码：`src/bwt.rs`
- 基准夹具：`/tmp/minibwa-ref.fa`（1 Mbp）、`/tmp/minibwa-bench-20000.fq`
- 本记录：`changelog/2026-11-08-bwt-rank-no-full-copy.md`
