#!/bin/bash
# 对比 C 版 minibwa 与 minibwa-rust 的 SAM 输出。
# 用法: ./scripts/compare.sh <reference.fa> <reads1.fq> [reads2.fq]
#
# 参考二进制默认 ../minibwa/minibwa（相对仓库根），可用 C_MINIBWA 覆盖。
# 两端都用当前 CLI：index <fa> <prefix> ； map [-t1] <prefix> <reads...>

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
C_BIN="${C_MINIBWA:-$REPO_ROOT/../minibwa/minibwa}"
RUST_BIN="${RUST_BIN:-$REPO_ROOT/target/release/minibwa-rust}"

REFERENCE="${1:?Usage: $0 <reference.fa> <reads1.fq> [reads2.fq]}"
READS1="${2:?Usage: $0 <reference.fa> <reads1.fq> [reads2.fq]}"
READS2="${3:-}"

if [[ ! -x "$C_BIN" ]]; then
    echo "ERROR: C minibwa 未找到或不可执行: $C_BIN" >&2
    echo "       请 clone https://github.com/lh3/minibwa 到 ../minibwa 并 make。" >&2
    exit 1
fi

WORKDIR="${COMPARE_WORKDIR:-$(pwd)}"
C_PREFIX="$WORKDIR/compare_c"
R_PREFIX="$WORKDIR/compare_r"
C_SAM="$WORKDIR/minibwa_c.sam"
R_SAM="$WORKDIR/minibwa_rust.sam"

echo "=== Build minibwa-rust (release) ==="
(cd "$REPO_ROOT" && cargo build --release)
if [[ ! -x "$RUST_BIN" ]]; then
    echo "ERROR: 未找到 $RUST_BIN" >&2
    exit 1
fi

echo "=== Index (C) ==="
"$C_BIN" index -t1 "$REFERENCE" "$C_PREFIX"

echo "=== Index (Rust) ==="
"$RUST_BIN" index -t1 "$REFERENCE" "$R_PREFIX"

echo "=== Map (C) ==="
if [[ -n "$READS2" ]]; then
    "$C_BIN" map -t1 "$C_PREFIX" "$READS1" "$READS2" > "$C_SAM"
else
    "$C_BIN" map -t1 "$C_PREFIX" "$READS1" > "$C_SAM"
fi

echo "=== Map (Rust) ==="
if [[ -n "$READS2" ]]; then
    "$RUST_BIN" map -t1 "$R_PREFIX" "$READS1" "$READS2" > "$R_SAM"
else
    "$RUST_BIN" map -t1 "$R_PREFIX" "$READS1" > "$R_SAM"
fi

count_mapped() {
    local sam_file="$1"
    local total unmapped mapped
    total=$(grep -cv '^@' "$sam_file" || true)
    unmapped=$(awk '!/^@/ && and($2, 4)' "$sam_file" | wc -l || true)
    mapped=$((total - unmapped))
    echo "  Total records: $total"
    echo "  Mapped:        $mapped"
    echo "  Unmapped:      $unmapped"
}

echo ""
echo "=== Comparison ==="
echo "--- C minibwa ($C_BIN) ---"
count_mapped "$C_SAM"
echo "--- minibwa-rust ---"
count_mapped "$R_SAM"

c_body=$(grep -cv '^@' "$C_SAM" || true)
r_body=$(grep -cv '^@' "$R_SAM" || true)
if [[ "$c_body" == "$r_body" ]]; then
    echo "Record count: match ($c_body)"
else
    echo "Record count: differ (C=$c_body Rust=$r_body)"
    exit 1
fi

echo ""
echo "SAM files: $C_SAM  $R_SAM"
echo "去 header 后可用: diff <(grep -v '^@' $C_SAM) <(grep -v '^@' $R_SAM)"
