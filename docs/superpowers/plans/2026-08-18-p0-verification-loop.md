# P0 验证闭环 Implementation Plan

> **For agentic workers:** 本会话按用户指示顺序执行（inline）。后续 P1/P2/P3 各自独立成计划。

**Goal:** 让 C/Rust 字节级对照在本机和 CI 上成为强制、可观察的事实，并修掉 `-x sr` 误置 ADAP。

**Architecture:** 将 C 参考 pin 到固定 commit；parity 测试在缺少 C 二进制时失败（CI）或仍可本地 skip；文档与 `scripts/compare.sh` 只描述真实接口。

**Tech Stack:** Rust 1.85, GitHub Actions, C minibwa (`make`), chrM 测试数据。

**Spec:** 对话中已确认的 P0 清单；`openspec/specs/reference-parity/spec.md`。

## Global Constraints

- 纯 Rust，`unsafe_code = "forbid"`
- 文档只保留中文
- 字段命名与 C 版一致
- 凡与 C 交互的格式/算法必须有对照测试
- 不猜测未 pin 的 C 行为

---

### Task 1: Pin C 参考源

**Files:**
- Create: `/home/shane/github/open-genomics/minibwa`（仓库外 clone，不进 git）
- Modify: `AGENTS.md`, `README.md`, `openspec/project.md`

- [ ] Clone `https://github.com/lh3/minibwa` 到 `../minibwa`，记录 HEAD SHA
- [ ] `make` 编译出 `../minibwa/minibwa`
- [ ] 文档写明 pin commit、构建命令、parity 前提

### Task 2: `-x sr` 不置 ADAP

**Files:**
- Modify: `src/opt.rs`
- Test: `src/opt.rs` 模块测试

- [ ] 失败测试：`apply_preset("sr")` 后 `flag & ADAP == 0`，`apply_preset("adap")` 后非 0
- [ ] 删掉 sr 分支里无条件的 `self.flag |= flag::ADAP`

### Task 3: 单端 SAM/PAF parity

**Files:**
- Modify: `tests/c_parity.rs`

- [ ] `se_sam_parity`：chrM read_1，去 header 后逐行
- [ ] `se_paf_parity`：`-f` PAF 逐字节

### Task 4: compare.sh

**Files:**
- Modify: `scripts/compare.sh`, `data/README.md`

- [ ] 调用 `../minibwa/minibwa` 与 `cargo run --release -- map`
- [ ] 索引 `.l2b`/`.mbw`；不再使用 `bwa` 或 `align -i .fm`

### Task 5: CI 强制 parity

**Files:**
- Modify: `.github/workflows/ci.yml`

- [ ] checkout pin 的 C commit，安装 zlib/make，`make`
- [ ] 将二进制放到 `../minibwa/minibwa`（相对 crate root）
- [ ] `cargo test` 必须真正执行 c_parity，缺二进制即失败

### Task 6: 文档同步

**Files:**
- Modify: `README.md`, `AGENTS.md`, `data/README.md`

- [ ] 里程碑表反映 P0 状态
- [ ] 示例命令与当前 CLI 一致
