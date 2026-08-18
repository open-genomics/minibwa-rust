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

- **C minibwa**: `https://github.com/lh3/minibwa`, commit in `C_MINIBWA_COMMIT`.
  Local path: `../minibwa`（二进制 `../minibwa/minibwa`）。CI 编译该 commit 并设置 `MINIBWA_REQUIRE_PARITY=1`。
- **`scripts/compare.sh`**: C `minibwa map` vs Rust `map`；索引 `.l2b`/`.mbw`。

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
