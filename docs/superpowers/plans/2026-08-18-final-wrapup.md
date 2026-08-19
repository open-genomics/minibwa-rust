# Final Wrap-up Plan

> **Target executor:** 便宜模型（如 Claude Haiku / GPT-4o-mini / 本地小模型）
> **原则：** 每步独立可验证，不跨步依赖；只做收尾不做新功能。

---

## 背景

minibwa-rust 是 C 版 minibwa 的 Rust 重写，核心算法（索引构建、SMEM 种子、链化、ksw2 比对、配对端、SAM/PAF 输出）已与 C 版逐字节/逐行对照通过。当前仓库状态：

- 64 个 Rust 单元测试 + 16 个 C parity 测试 + 1 个 doc-test，全部通过
- CI 已配置 C 参考二进制构建与 parity 强制检查
- Release workflow 已配置多平台构建
- `scripts/compare.sh` 已在最新 commit 中修复为使用 `minibwa`（而非 `bwa`），但缺少退出码

**不在本计划范围内的已知未做项：**
- SIMD ksw2（受 `unsafe_code = "forbid"` 阻塞）
- 读/比对/写流水线（纯性能优化，非正确性）
- `map --meth` BS-seq 映射（`index --meth` 已完成，映射未做）

这些已在 README 中标注为"未做"，本计划不再处理。

---

## Task 1: 给 `scripts/compare.sh` 添加退出码

**目标：** 脚本当前已使用正确的 CLI（`minibwa` + `map`），但 record count 不一致时不会以非零退出码退出，导致 CI 或自动化脚本无法检测失败。

**文件：** `scripts/compare.sh`

**具体改动：**

找到第 77-81 行：
```bash
if [[ "$c_body" == "$r_body" ]]; then
    echo "Record count: match ($c_body)"
else
    echo "Record count: differ (C=$c_body Rust=$r_body)"
fi
```

替换为：
```bash
if [[ "$c_body" == "$r_body" ]]; then
    echo "Record count: match ($c_body)"
else
    echo "Record count: differ (C=$c_body Rust=$r_body)"
    exit 1
fi
```

**验证：**
```bash
cd /home/shane/github/open-genomics/minibwa-rust
cargo build --release
./scripts/compare.sh data/toy.fa data/toy_reads.fq
echo $?  # 应输出 0
```

---

## Task 2: 归档已完成的 change proposal

**目标：** `openspec/changes/correct-parity-contract-docs/` 已标记为 `Completed`，应移入 archive。

**文件操作：**

```bash
cd /home/shane/github/open-genomics/minibwa-rust
mv openspec/changes/correct-parity-contract-docs openspec/changes/archive/
```

然后更新 `openspec/changes/archive/` 目录下的 `correct-parity-contract-docs/proposal.md`：
- 将第 5 行 `Status: Applying` 改为 `Status: Archived`

**验证：**
```bash
ls openspec/changes/archive/correct-parity-contract-docs/
# 应看到 proposal.md, design.md, tasks.md, verification.md, specs/
```

---

## Task 3: 更新 README 最终状态

**目标：** 让 README 准确反映项目最终交付状态，明确区分"已完成"、"已知限制"、"不再计划"。

**文件：** `README.md`

**具体改动：**

1. 将第 34-35 行（未做项）扩展为更清晰的"已知限制"小节：

```markdown
## 已知限制

以下功能在 C 版 minibwa 中存在，但本 Rust 重写因政策或优先级原因不实现：

| 功能 | 原因 |
|------|------|
| SIMD ksw2（SSE/NEON 条带 DP） | 与 `unsafe_code = "forbid"` 冲突；标量版已与 C 对照一致 |
| 读/比对/写流水线 | 纯性能优化，不影响正确性；当前 batch 内 rayon 并行已足够 |
| `map --meth` BS-seq 映射 | `index --meth` 已与 C 字节一致，映射路径未实现 |
```

2. 在"当前状态"表格后添加一行：

```markdown
| — | 最终收尾：compare.sh 修复、文档定型 | ✅ |
```

3. 将第 7 行改为：

```markdown
- 索引文件与 C 版 minibwa **字节级兼容**（`.l2b` 参考序列 + `.mbw` BWT/采样 SA；`--meth` 时另有 `.meth.mbw`）
- 算法逐模块对照移植，以 **已 pin 的 C 版 commit** 输出为基准做一致性验证
- 纯 Rust 实现（后缀数组构建使用 libsais 的纯 Rust 移植 `libsais-rs`），无 C 依赖
- **本项目已功能完备（标准 DNA 比对），性能优化与 BS-seq 映射不在范围内**
```

**验证：**
```bash
cat README.md | grep -A5 "已知限制"
cat README.md | grep "功能完备"
```

---

## Task 4: 最终验证

**目标：** 确保所有改动后项目仍通过全部检查。

**命令：**
```bash
cd /home/shane/github/open-genomics/minibwa-rust

# 1. 格式化检查
cargo fmt --all -- --check

# 2. Clippy
cargo clippy --all-targets --all-features -- -D warnings

# 3. 测试（含 C parity）
cargo test --all-targets --all-features

# 4. 确认 compare.sh 可运行
./scripts/compare.sh data/toy.fa data/toy_reads.fq
```

**预期结果：** 全部通过，exit code 0。

---

## 执行顺序

1. Task 1（修复 compare.sh）
2. Task 2（归档 change proposal）
3. Task 3（更新 README）
4. Task 4（最终验证）

每步完成后执行该步的验证命令，确认通过后再继续。

---

## 不做的内容（明确排除）

- ❌ 不实现 SIMD ksw2
- ❌ 不实现读/比对/写流水线
- ❌ 不实现 `map --meth`
- ❌ 不修改 Rust 算法代码
- ❌ 不添加新的 parity 测试
- ❌ 不修改 CI/CD 配置
- ❌ 不创建 CHANGELOG
- ❌ 不提交代码（由用户自行 git commit）