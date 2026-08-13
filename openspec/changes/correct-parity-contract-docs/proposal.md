# Change Proposal: correct-parity-contract-docs

## Metadata

- Status: `Applying`
- Repository: `open-genomics/minibwa-rust`
- Audit base: `dad6e99b064a93bd017ad819b529686b46413b02`
- Capability: `reference-parity`
- Task IDs: `MBW-DOC-001`, `ORG-GOV-001`
- Decision IDs: none

## Why

核心承诺是与 C minibwa 一致，但 `tests/c_parity.rs` 在参考二进制或测试数据缺失时返回成功，当前 CI 未准备这些依赖；README 仍容易让用户把 parity 用例存在误解为 CI 已验证。`scripts/compare.sh` 还使用漂移的命令/文件假设。

## Changes

**Parity claim visibility**

- From: 参考依赖缺失会静默 skip，文档没有清楚区分“测试代码存在”和“本次运行实际验证”。
- To: 文档准确说明 reference 条件、当前 CI 证明范围和 skip；compare 脚本只保留已验证接口，否则明确退役/阻塞。
- Reason: 不把未执行的核心测试表述为兼容证明。
- Impact: 文档/工具真实性修复；不固定未知 C commit，不改变 Rust 算法。

## Scope

- `README.md`、`AGENTS.md`、里程碑和 parity 说明；
- `scripts/compare.sh` 的当前接口核查与最小修复，若不依赖未知来源；
- 仓库内 `reference-parity` spec；
- 允许为 skip 条件添加可观察输出/测试，但不把 CI 改成强制模式。

## Out of scope

- 不改变已声明的 canonical C URL，也不猜测参考 commit；
- 不自动下载 upstream；
- 不实施 `MBW-REF-001/CI-001/FIXTURE-001/REL-001`；
- 不接受或更新 parity output；
- 不修改 alignment/index 算法。

## Blocker rule

如果修复 compare script 必须知道尚未确认的 C commit/build/interface，则保留脚本修改 task 为 Blocked，只完成准确文档，不用 `bwa` 或假数据替代 minibwa。

## Approval

- Apply approval: `authorized by organization owner`
