# Tasks: correct-parity-contract-docs

## 1. Baseline

- [x] 1.1 记录 HEAD/status，运行 cargo test 并确认 live parity 是否执行
  - HEAD: dad6e99, clean; C binary present, 8 parity tests executed and passed
- [x] 1.2 核对 tests/c_parity.rs 的全部 skip 条件和 CI 提供的依赖
  - Skip condition: `../minibwa/minibwa` not found; CI does not build C binary
- [x] 1.3 逐条核对 scripts/compare.sh 使用的命令、选项和产物是否存在
  - Script uses `bwa` (not `minibwa`); Rust CLI uses `map` not `align`; no `-o` flag; index ext `.l2b`/`.mbw` not `.fm`

## 2. Truthful contract

- [x] 2.1 建立 openspec/project.md 与 reference-parity change
- [x] 2.2 更新 README/AGENTS 区分 Rust tests、live parity、CI scope
  - README: parity prerequisites note; AGENTS.md: verification section updated
- [x] 2.3 让 skip 原因可观察
  - tests/c_parity.rs: skip message updated to `[PARITY SKIP]` prefix
- [x] 2.4 compare script 接口漂移已记录，修复 Blocked（需确认 C source）
  - scripts/compare.sh: interface drift documented in header and runtime warning
- [x] 2.5 不存在 CHANGELOG，无需更新

## 3. Verification

- [x] 3.1 cargo fmt --all -- --check - passed
- [x] 3.2 cargo clippy --all-targets --all-features -- -D warnings - passed
- [x] 3.3 cargo test --all-targets --all-features - 70 tests passed (62 unit + 8 parity)
- [x] 3.4 有参考依赖时 parity 执行通过；无依赖时 skip 消息可观察
- [x] 3.5 git diff --check, scope audit, verification.md - passed
