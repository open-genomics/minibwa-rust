# 变更记录：seed_intv 第二轮 SMEM 以父区间右端为界

- 日期：2026-11-08
- 变更类型：算法修复（与 C 版 seed.c 语义对齐）

## 内容
- C 版 `seed.c:47` 第二轮调用 `mb_bwt_smem(bwt, en, seq, x, sub_min_len, ...)`，
  以父区间右端 `en` 作为搜索长度；`bwt.c` 的 `mb_bwt_smem` 在 `len - x < min_len`
  时直接返回 `len`，前向扩展循环上界也是 `len`。
- Rust 的 `Bwt::smem` 内部固定使用 `q.len()` 作为长度，而 `seed_intv` 第二轮
  原先传入整条 `seq`，导致第二轮搜索越过父区间右端，产生 C 版不存在的子区间。
- 修复：第二轮把查询截断到父区间 `&seq[..en]` 再传给 `smem`，与 C 的 `len=en`
  完全等价（`en` 由第一轮 SMEM 产生，恒 ≤ `seq.len()`）。

## 影响
- 含重复区的 reads 不再产生多余子区间；这些区间此前会进入 anchors，进而可能
  改变链化/hit 结果。
- 第一轮 SMEM 与不满足第二轮条件的父区间行为不变。
- 未传 `-Y` 等其余路径无变化。

## 验证
- 新增 `src/seed.rs::tests::seed_intv_pass2_uses_parent_end`：24bp 参考
  `CAAAAACAAAAACCAACAACAAAA`、query `CCAACAACAAAAA`（nt4）、min_len=6、
  max_sub_occ=10。修复前出现 `info=(7<<32)|13, size=2` 的区间，修复后消失。
- 仓库外临时校验程序（/tmp/minicheck.rs，未入库）按 C 语义实现 seed_intv：
  修复后 Rust 输出与 C 语义逐项相等（equal = true）；修复前多出该区间。
- `cargo test -j 4`：88 passed / 0 failed（65 lib + 1 bin + 16 c_parity + 5
  test_corrupt_index + 1 doctest）；`cargo clippy --all-targets -- -D warnings`、
  `cargo fmt --all -- --check` 通过。
- 本地缺 `../minibwa` 与 chrM 数据，16 个 c_parity 测试全部打印
  `[PARITY SKIP]`，本次未验证与 C 二进制的字节级输出一致性。

## 路径
- 代码：`src/seed.rs`
- 测试：`src/seed.rs`（tests 模块）
- 本记录：`changelog/2026-11-08-seed-pass2-parent-end.md`
