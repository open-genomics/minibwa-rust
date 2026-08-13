# minibwa-rust - Project Context

## Identity

- **Canonical repository**: `open-genomics/minibwa-rust`
- **Lifecycle**: experimental; Rust rewrite of C minibwa
- **Language**: Rust (unsafe_code = "forbid")

## Core contracts

| Capability | Path | Description |
|---|---|---|
| `reference-parity` | `openspec/specs/reference-parity/` | C minibwa parity claims, skip behavior, comparison script |

## External boundaries

- **C minibwa**: reference implementation at `https://github.com/lh3/minibwa`.
  Local path: `../minibwa/minibwa`. No pinned commit confirmed.
- **`bwa` (Heng Li)**: used by `scripts/compare.sh`, not `minibwa`. Script has interface drift.

## Validation commands

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
```

## Authority rules

- `src/` is the implementation source of truth.
- Models must not commit, push, create PRs, or publish without explicit authorization.
- Parity claims must distinguish Rust-only tests from C reference comparisons.
