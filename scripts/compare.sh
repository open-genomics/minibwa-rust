#!/bin/bash
# 对比原版 BWA 与 minibwa-rust 的对齐结果
# 用法: ./scripts/compare.sh <reference.fa> <reads.fq>
#
# 前置条件: 系统上安装了原版 bwa (bwa mem)
# 输出: 在当前目录生成 bwa_result.sam 和 minibwa_rust_result.sam，并打印简要统计对比

set -euo pipefail

REFERENCE="${1:?Usage: $0 <reference.fa> <reads.fq>}"
READS="${2:?Usage: $0 <reference.fa> <reads.fq>}"

BWA_SAM="bwa_result.sam"
RUST_SAM="minibwa_rust_result.sam"
INDEX_PREFIX="compare_ref"

echo "=== Step 1: Build BWA index ==="
if command -v bwa &>/dev/null; then
    bwa index "$REFERENCE"
else
    echo "WARNING: bwa not found in PATH, skipping C BWA comparison"
fi

echo "=== Step 2: Build minibwa-rust index ==="
cargo run --release -- index "$REFERENCE" -o "$INDEX_PREFIX"

echo "=== Step 3: Run BWA MEM ==="
if command -v bwa &>/dev/null; then
    bwa mem "$REFERENCE" "$READS" > "$BWA_SAM" 2>/dev/null
    echo "BWA result: $BWA_SAM"
else
    echo "Skipped (bwa not found)"
fi

echo "=== Step 4: Run minibwa-rust align ==="
cargo run --release -- align -i "${INDEX_PREFIX}.fm" "$READS" -o "$RUST_SAM"
echo "minibwa-rust result: $RUST_SAM"

echo ""
echo "=== Comparison ==="

count_mapped() {
    local sam_file="$1"
    # 跳过 header 行，统计 mapped (FLAG 没有 0x4 位) 和 unmapped
    local total mapped unmapped
    total=$(grep -cv '^@' "$sam_file" || true)
    unmapped=$(awk '!/^@/ && and($2, 4)' "$sam_file" | wc -l || true)
    mapped=$((total - unmapped))
    echo "  Total reads: $total"
    echo "  Mapped:      $mapped"
    echo "  Unmapped:    $unmapped"
}

if [ -f "$BWA_SAM" ]; then
    echo "--- BWA MEM ---"
    count_mapped "$BWA_SAM"
fi

echo "--- minibwa-rust ---"
count_mapped "$RUST_SAM"

echo ""
echo "Done. You can manually diff the SAM files for detailed comparison."
