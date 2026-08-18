# 测试数据与基准结果

## 仓库内数据

- `toy.fa` — 小型参考序列，包含 3 条 contig（chr1, chr2, chr3）
- `toy_reads.fq` — 对应 reads，覆盖精确匹配、错配、indel、反向互补、unmapped 等情况

chrM 对照数据在 C 仓库 `../minibwa/test/`（`chrM-human.fa.gz`、`chrM-read_1.fa.gz`、`chrM-read_2.fa.gz`），不复制进本仓库。

## 与 C 版 minibwa 对照

```bash
# 1. 准备 C 二进制（commit 见仓库根 C_MINIBWA_COMMIT）
git clone https://github.com/lh3/minibwa.git ../minibwa
make -C ../minibwa

# 2. 建索引并对齐
cargo run --release -- index data/toy.fa data/toy
cargo run --release -- map data/toy data/toy_reads.fq > data/toy_rust.sam

# 3. 一键对比（C map vs Rust map）
./scripts/compare.sh data/toy.fa data/toy_reads.fq
```
