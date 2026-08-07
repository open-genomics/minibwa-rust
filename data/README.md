# 测试数据与基准结果

## 参考数据

- `toy.fa` — 小型参考序列，包含 3 个 contig（chr1, chr2, chr3）
- `toy_reads.fq` — 对应 reads，覆盖精确匹配、错配、indel、反向互补、unmapped 等情况

## 使用 C 版 BWA 生成基准结果

```bash
# 1. 构建索引
bwa index toy.fa

# 2. 运行 BWA-MEM
bwa mem toy.fa toy_reads.fq > toy_bwa_mem.sam

# 3. 查看结果
cat toy_bwa_mem.sam
```

## 使用 minibwa-rust 生成结果

```bash
# 1. 构建索引
cargo run --release -- index data/toy.fa -o data/toy

# 2. 运行对齐
cargo run --release -- align -i data/toy.fm data/toy_reads.fq -o data/toy_rust.sam
```

## 对比方法

使用 `scripts/compare.sh` 脚本可以快速对比两者的 mapped/unmapped 数量：

```bash
./scripts/compare.sh data/toy.fa data/toy_reads.fq
```
