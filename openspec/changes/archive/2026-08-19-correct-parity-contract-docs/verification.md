# Verification: correct-parity-contract-docs

- Status: `Completed`
- Ready to archive: `yes`
- Verifier: `implementing agent (self-verified)`
- Date: 2026-08-13

## Environment

- HEAD: `dad6e99b064a93bd017ad819b529686b46413b02` (matches audit base)
- Working tree: clean before apply
- C minibwa binary: present at `../minibwa/minibwa` (parity tests executed, not skipped)

## Requirement -> Evidence Matrix

| Requirement | Scenario | Evidence | Result |
|---|---|---|---|
| Parity claims identify reference prerequisites | User reads validation instructions | README: parity prerequisites note added; AGENTS.md: verification section updated | passed |
| Current CI scope is truthful | CI runs without C reference | README: "当前 CI 未构建 C 参考二进制，CI 通过仅表示 Rust 自身测试通过" | passed |
| Skip behavior is observable | Reference binary missing | `tests/c_parity.rs`: skip message updated to `[PARITY SKIP]` prefix with clear explanation | passed |
| Comparison helper matches supported interfaces | Helper cannot be validated | `scripts/compare.sh`: interface drift documented in header and runtime warning; uses `bwa` not `minibwa`; Rust CLI commands outdated | passed (documented as unsupported) |

## Command Results

| Command | Exit status | Summary |
|---|---|---|
| `cargo fmt --all -- --check` | 0 | No formatting issues |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | No warnings |
| `cargo test --all-targets --all-features` | 0 | 62 unit + 8 parity tests passed (C binary present) |
| `git diff --check` | 0 | No whitespace errors |

## Notes

- C minibwa binary exists in this workspace, so parity tests executed and passed.
- CI does not build C binary; documentation now truthfully states this.
- `scripts/compare.sh` uses `bwa` (not `minibwa`) and has outdated Rust CLI commands; documented as having interface drift, not fixed (blocked on confirming C source).
