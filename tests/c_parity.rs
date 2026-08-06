//! 与 C 版 minibwa 的字节级对照测试。
//!
//! 需要相邻目录存在已编译的 C 版 minibwa（`../minibwa/minibwa`）；
//! 不存在时跳过并打印提示（CI 无 C 二进制时不视为失败）。

use std::path::{Path, PathBuf};
use std::process::Command;

const TOY_FA: &str = r#">chr1 human mitochondrial-like test
ACGTACGTACGTNNNNacgtACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2  second contig with Ns and lowercase
TTTTNNNNGGGGCCCCaaaaccccAAAACCCCggggnnnnTTTTAAAA
>chr3
ACACACACACACACACACACGTGTGTGTGTGTGTGTGTGTGTGTGTGT
"#;

fn c_minibwa() -> Option<PathBuf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../minibwa/minibwa");
    if p.exists() {
        Some(p)
    } else {
        eprintln!("跳过 C 版对照测试：未找到 {}", p.display());
        None
    }
}

fn rust_bin() -> &'static str {
    env!("CARGO_BIN_EXE_minibwa-rust")
}

fn write(path: &Path, content: &str) {
    std::fs::write(path, content).unwrap();
}

fn run(cmd: &mut Command) {
    let st = cmd.status().unwrap();
    assert!(st.success(), "命令失败: {cmd:?}");
}

/// 用 C/Rust 各自建索引并逐字节比较 .l2b 与 .mbw。
fn parity_index(dir: &Path, c: &Path, extra_args: &[&str]) {
    let fa = dir.join("toy.fa");
    write(&fa, TOY_FA);

    let mut cc = Command::new(c);
    cc.arg("index").args(extra_args).arg("-t1").arg(&fa).arg(dir.join("c"));
    run(&mut cc);

    let mut rc = Command::new(rust_bin());
    rc.arg("index")
        .args(extra_args)
        .arg("-t1")
        .arg(&fa)
        .arg(dir.join("r"));
    run(&mut rc);

    let c_l2b = std::fs::read(dir.join("c.l2b")).unwrap();
    let r_l2b = std::fs::read(dir.join("r.l2b")).unwrap();
    assert_eq!(c_l2b, r_l2b, ".l2b 字节不一致 (args={extra_args:?})");

    let c_mbw = std::fs::read(dir.join("c.mbw")).unwrap();
    let r_mbw = std::fs::read(dir.join("r.mbw")).unwrap();
    assert_eq!(c_mbw, r_mbw, ".mbw 字节不一致 (args={extra_args:?})");
}

#[test]
fn index_parity_default() {
    let Some(c) = c_minibwa() else { return };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-default");
    std::fs::create_dir_all(&dir).unwrap();
    parity_index(&dir, &c, &[]);
}

#[test]
fn index_parity_sa_bit_and_seed() {
    let Some(c) = c_minibwa() else { return };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-u2s42");
    std::fs::create_dir_all(&dir).unwrap();
    parity_index(&dir, &c, &["-u", "2", "-s", "42"]);
}

#[test]
fn getref_parity() {
    let Some(c) = c_minibwa() else { return };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-getref");
    std::fs::create_dir_all(&dir).unwrap();
    parity_index(&dir, &c, &[]);

    // 交叉验证：Rust getref 读 C 建的 .l2b
    let out = Command::new(rust_bin())
        .arg("getref")
        .arg(dir.join("c.l2b"))
        .output()
        .unwrap();
    assert!(out.status.success());

    let c_out = Command::new(&c)
        .arg("getref")
        .arg(dir.join("c.l2b"))
        .output()
        .unwrap();
    assert!(c_out.status.success());
    assert_eq!(c_out.stdout, out.stdout, "getref 输出不一致");
}

#[test]
fn rust_loads_c_index_chrm() {
    // 加载 C 版对真实数据（chrM）建的索引：验证 gzip FASTA 路径与大一点的基因组
    let Some(c) = c_minibwa() else { return };
    let mb = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../minibwa/test/chrM-human.fa.gz");
    if !mb.exists() {
        eprintln!("跳过：缺少 {}", mb.display());
        return;
    }
    let dir = std::env::temp_dir().join("minibwa-rust-parity-chrm");
    std::fs::create_dir_all(&dir).unwrap();

    let mut cc = Command::new(c);
    cc.arg("index").arg("-t1").arg(&mb).arg(dir.join("c"));
    run(&mut cc);
    let mut rc = Command::new(rust_bin());
    rc.arg("index").arg("-t1").arg(&mb).arg(dir.join("r"));
    run(&mut rc);

    assert_eq!(
        std::fs::read(dir.join("c.l2b")).unwrap(),
        std::fs::read(dir.join("r.l2b")).unwrap(),
        "chrM .l2b 不一致"
    );
    assert_eq!(
        std::fs::read(dir.join("c.mbw")).unwrap(),
        std::fs::read(dir.join("r.mbw")).unwrap(),
        "chrM .mbw 不一致"
    );
}
