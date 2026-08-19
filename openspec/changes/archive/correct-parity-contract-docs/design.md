# Design: correct-parity-contract-docs

## Evidence

- `tests/c_parity.rs` 检查相邻 C binary/数据，缺失时测试成功退出；
- `.github/workflows/ci.yml` 运行 Cargo tests，但未准备 C 参考构建和 chrM 输入；
- release workflow 同样会经过可 skip 路径；
- `scripts/compare.sh` 与当前 Rust CLI/产物存在漂移。

## Target support statement

明确区分：

1. Rust-only tests：CI 必然运行；
2. live C parity tests：只有 reference binary/data 可用时运行，当前 CI 不构成证明；
3. frozen reference fixtures：尚未建立；
4. release parity gate：尚未建立。

## Allowed surface

- `openspec/`
- `README.md`, `AGENTS.md`
- `scripts/compare.sh` 及其 shell test（仅当可由当前仓库事实确定）
- tests 中只与 skip 可观察性相关的最小更改
- `CHANGELOG.md`

## Follow-up boundary

未来强制 parity 使用仓库已声明的 `https://github.com/lh3/minibwa`，但必须先确认完整 commit、build command 和 fixture provenance，再进入 `pin-c-reference-source`。本 change 不降低这个前置条件。
