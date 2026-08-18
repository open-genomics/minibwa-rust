//! 与 C 版 minibwa 的字节级对照测试。
//!
//! 需要相邻目录存在已编译的 C 版 minibwa（`../minibwa/minibwa`），commit 见仓库根目录
//! `C_MINIBWA_COMMIT`。默认缺失时跳过并打印 `[PARITY SKIP]`；设置
//! `MINIBWA_REQUIRE_PARITY=1`（CI）时缺失即失败。

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
        return Some(p);
    }
    let msg = format!(
        "C 版 minibwa 二进制未找到 ({})。Parity 测试未执行，本次运行不验证 C 兼容性。",
        p.display()
    );
    if std::env::var_os("MINIBWA_REQUIRE_PARITY").is_some() {
        panic!("MINIBWA_REQUIRE_PARITY=1：{msg}");
    }
    eprintln!("\n[PARITY SKIP] {msg}\n");
    None
}

fn strip_sam_header(out: &[u8]) -> Vec<String> {
    std::str::from_utf8(out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(ToString::to_string)
        .collect()
}

fn require_chrm(name: &str) -> Option<PathBuf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../minibwa/test").join(name);
    if p.exists() {
        Some(p)
    } else if std::env::var_os("MINIBWA_REQUIRE_PARITY").is_some() {
        panic!("MINIBWA_REQUIRE_PARITY=1：缺少 {}", p.display());
    } else {
        eprintln!("跳过：缺少 {}", p.display());
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
    rc.arg("index").args(extra_args).arg("-t1").arg(&fa).arg(dir.join("r"));
    run(&mut rc);

    let c_l2b = std::fs::read(dir.join("c.l2b")).unwrap();
    let r_l2b = std::fs::read(dir.join("r.l2b")).unwrap();
    assert_eq!(c_l2b, r_l2b, ".l2b 字节不一致 (args={extra_args:?})");

    let c_mbw = std::fs::read(dir.join("c.mbw")).unwrap();
    let r_mbw = std::fs::read(dir.join("r.mbw")).unwrap();
    assert_eq!(c_mbw, r_mbw, ".mbw 字节不一致 (args={extra_args:?})");

    if extra_args.contains(&"--meth") {
        let c_meth = std::fs::read(dir.join("c.meth.mbw")).unwrap();
        let r_meth = std::fs::read(dir.join("r.meth.mbw")).unwrap();
        assert_eq!(c_meth, r_meth, ".meth.mbw 字节不一致 (args={extra_args:?})");
    }
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
fn index_parity_meth() {
    let Some(c) = c_minibwa() else { return };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-meth");
    std::fs::create_dir_all(&dir).unwrap();
    parity_index(&dir, &c, &["--meth"]);
}

#[test]
fn getref_parity() {
    let Some(c) = c_minibwa() else { return };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-getref");
    std::fs::create_dir_all(&dir).unwrap();
    parity_index(&dir, &c, &[]);

    // 交叉验证：Rust getref 读 C 建的 .l2b
    let out = Command::new(rust_bin()).arg("getref").arg(dir.join("c.l2b")).output().unwrap();
    assert!(out.status.success());

    let c_out = Command::new(&c).arg("getref").arg(dir.join("c.l2b")).output().unwrap();
    assert!(c_out.status.success());
    assert_eq!(c_out.stdout, out.stdout, "getref 输出不一致");
}

#[test]
fn rust_loads_c_index_chrm() {
    // 加载 C 版对真实数据（chrM）建的索引：验证 gzip FASTA 路径与大一点的基因组
    let Some(c) = c_minibwa() else { return };
    let Some(mb) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
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

#[test]
fn fastmap_parity() {
    // 种子级对照：C/Rust 的 fastmap 输出必须逐字节一致
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-fastmap");
    std::fs::create_dir_all(&dir).unwrap();

    run(Command::new(&c).arg("index").arg("-t1").arg(&fa).arg(dir.join("c")));
    run(Command::new(rust_bin()).arg("index").arg("-t1").arg(&fa).arg(dir.join("r")));

    for args in [
        &["-l", "19", "-s", "1", "-w", "20"][..],
        &["-l", "15", "-s", "2", "-w", "5"][..],
    ] {
        let c_out = Command::new(&c)
            .arg("fastmap")
            .args(args)
            .arg(dir.join("c"))
            .arg(&reads)
            .output()
            .unwrap();
        let r_out = Command::new(rust_bin())
            .arg("fastmap")
            .args(args)
            .arg(dir.join("r"))
            .arg(&reads)
            .output()
            .unwrap();
        assert!(c_out.status.success() && r_out.status.success());
        assert_eq!(c_out.stdout, r_out.stdout, "fastmap 输出不一致 (args={args:?})");
    }
}

#[test]
fn chain_only_paf_parity() {
    // 链级对照：--chain-only PAF 输出必须逐字节一致（覆盖链化/hit/mapQ）
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_2.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-chain");
    std::fs::create_dir_all(&dir).unwrap();

    run(Command::new(&c).arg("index").arg("-t1").arg(&fa).arg(dir.join("c")));
    run(Command::new(rust_bin()).arg("index").arg("-t1").arg(&fa).arg(dir.join("r")));

    let cases: &[&[&str]] = &[
        &[],
        &["-x", "lr"],
        &["-x", "sr"],
        &["-k", "15", "-w", "50"],
        &["-m", "10", "-p", "0.3", "-N", "10"],
        &["-g", "500", "-W", "5000", "--outs", "0.5"],
        &["-x", "lr", "--outn", "3"],
        &["-m", "10", "-p", "0.3", "-N", "10", "-5"],
    ];
    for args in cases {
        let c_out = Command::new(&c)
            .arg("map")
            .arg("--chain-only")
            .arg("-t1")
            .args(*args)
            .arg(dir.join("c"))
            .arg(&reads)
            .output()
            .unwrap();
        let r_out = Command::new(rust_bin())
            .arg("map")
            .arg("--chain-only")
            .arg("-t1")
            .args(*args)
            .arg(dir.join("r"))
            .arg(&reads)
            .output()
            .unwrap();
        assert!(c_out.status.success() && r_out.status.success(), "args={args:?}");
        assert_eq!(c_out.stdout, r_out.stdout, "chain-only PAF 不一致 (args={args:?})");
    }
}

/// 单端 SAM：chrM read_1 去 header 后与 C 版逐行一致。
#[test]
fn se_sam_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-se-sam");
    std::fs::create_dir_all(&dir).unwrap();

    run(Command::new(&c).arg("index").arg("-t1").arg(&fa).arg(dir.join("c")));
    run(Command::new(rust_bin()).arg("index").arg("-t1").arg(&fa).arg(dir.join("r")));

    for args in [&["-t1"][..], &["-t1", "-x", "sr"][..], &["-t1", "-x", "lr"][..]] {
        let c_out = Command::new(&c)
            .arg("map")
            .args(args)
            .arg(dir.join("c"))
            .arg(&reads)
            .output()
            .unwrap();
        let r_out = Command::new(rust_bin())
            .arg("map")
            .args(args)
            .arg(dir.join("r"))
            .arg(&reads)
            .output()
            .unwrap();
        assert!(c_out.status.success() && r_out.status.success(), "args={args:?}");
        assert_eq!(
            strip_sam_header(&c_out.stdout),
            strip_sam_header(&r_out.stdout),
            "SE SAM 不一致 (args={args:?})"
        );
    }
}

/// 单端 PAF：`-f` 输出与 C 版逐字节一致。
#[test]
fn se_paf_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-se-paf");
    std::fs::create_dir_all(&dir).unwrap();

    run(Command::new(&c).arg("index").arg("-t1").arg(&fa).arg(dir.join("c")));
    run(Command::new(rust_bin()).arg("index").arg("-t1").arg(&fa).arg(dir.join("r")));

    let c_out = Command::new(&c)
        .arg("map")
        .arg("-f")
        .arg("-t1")
        .arg(dir.join("c"))
        .arg(&reads)
        .output()
        .unwrap();
    let r_out = Command::new(rust_bin())
        .arg("map")
        .arg("-f")
        .arg("-t1")
        .arg(dir.join("r"))
        .arg(&reads)
        .output()
        .unwrap();
    assert!(c_out.status.success() && r_out.status.success());
    assert_eq!(c_out.stdout, r_out.stdout, "SE PAF 不一致");
}

/// 配对端对照：PE SAM 输出必须与 C 版逐行一致（跳过 header）。
/// 覆盖：pestat 数据估计 → mb_pair → rescue → 配对 SAM（flag/RNEXT/PNEXT/TLEN/MC/MQ）。
#[test]
fn pe_sam_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(r1) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let Some(r2) = require_chrm("chrM-read_2.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-pe");
    std::fs::create_dir_all(&dir).unwrap();

    run(Command::new(&c).arg("index").arg("-t1").arg(&fa).arg(dir.join("c")));
    run(Command::new(rust_bin()).arg("index").arg("-t1").arg(&fa).arg(dir.join("r")));

    // 去 header 后逐行对比（@PG 的 CL 等元数据不同）
    for args in [&["-t1"][..], &["-t1", "-P"][..], &["-t1", "-5"][..]] {
        let c_out = Command::new(&c)
            .arg("map")
            .args(args)
            .arg(dir.join("c"))
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        let r_out = Command::new(rust_bin())
            .arg("map")
            .args(args)
            .arg(dir.join("r"))
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        assert!(c_out.status.success() && r_out.status.success(), "args={args:?}");
        assert_eq!(
            strip_sam_header(&c_out.stdout),
            strip_sam_header(&r_out.stdout),
            "PE SAM 不一致 (args={args:?})"
        );
    }
}

/// mem 命令配对端对照（bwa-mem 兼容 CLI，含 MD:Z）。
#[test]
fn mem_pe_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(r1) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let Some(r2) = require_chrm("chrM-read_2.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-mem");
    std::fs::create_dir_all(&dir).unwrap();

    run(Command::new(&c).arg("index").arg("-t1").arg(&fa).arg(dir.join("c")));
    run(Command::new(rust_bin()).arg("index").arg("-t1").arg(&fa).arg(dir.join("r")));

    for args in [&["-t1"][..], &["-t1", "-I", "300,80,500,100"][..], &["-t1", "-S"][..]] {
        let c_out = Command::new(&c)
            .arg("mem")
            .args(args)
            .arg(dir.join("c"))
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        let r_out = Command::new(rust_bin())
            .arg("mem")
            .args(args)
            .arg(dir.join("r"))
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap();
        assert!(c_out.status.success() && r_out.status.success(), "args={args:?}");
        assert_eq!(
            strip_sam_header(&c_out.stdout),
            strip_sam_header(&r_out.stdout),
            "mem PE SAM 不一致 (args={args:?})"
        );
    }
}

fn index_chrm(c: &Path, dir: &Path, fa: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    run(Command::new(c).arg("index").arg("-t1").arg(fa).arg(dir.join("c")));
    run(Command::new(rust_bin()).arg("index").arg("-t1").arg(fa).arg(dir.join("r")));
}

fn map_stdout(bin: &str, extra: &[&str], prefix: &Path, reads: &[&Path]) -> Vec<u8> {
    let mut cmd = Command::new(bin);
    cmd.arg("map").args(extra).arg("-t1").arg(prefix);
    for r in reads {
        cmd.arg(r);
    }
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "map 失败 extra={extra:?} stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

/// `map --eqx`：CIGAR 使用 =/X，与 C 版去 header 后逐行一致。
#[test]
fn map_eqx_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-eqx");
    index_chrm(&c, &dir, &fa);
    let extra = ["--eqx"];
    assert_eq!(
        strip_sam_header(&map_stdout(c.to_str().unwrap(), &extra, &dir.join("c"), &[&reads])),
        strip_sam_header(&map_stdout(rust_bin(), &extra, &dir.join("r"), &[&reads])),
        "map --eqx SAM 不一致"
    );
}

/// `map -b cs`：cs:Z 标签与 C 版逐行一致。
#[test]
fn map_cs_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-cs");
    index_chrm(&c, &dir, &fa);
    let extra = ["-b", "cs"];
    assert_eq!(
        strip_sam_header(&map_stdout(c.to_str().unwrap(), &extra, &dir.join("c"), &[&reads])),
        strip_sam_header(&map_stdout(rust_bin(), &extra, &dir.join("r"), &[&reads])),
        "map -b cs SAM 不一致"
    );
}

/// `map -b ds`：ds:Z 标签与 C 版逐行一致。
#[test]
fn map_ds_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-ds");
    index_chrm(&c, &dir, &fa);
    let extra = ["-b", "ds"];
    assert_eq!(
        strip_sam_header(&map_stdout(c.to_str().unwrap(), &extra, &dir.join("c"), &[&reads])),
        strip_sam_header(&map_stdout(rust_bin(), &extra, &dir.join("r"), &[&reads])),
        "map -b ds SAM 不一致"
    );
}

/// `map --hic` 等价于 `-5P`。
#[test]
fn map_hic_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(r1) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let Some(r2) = require_chrm("chrM-read_2.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-hic");
    index_chrm(&c, &dir, &fa);
    let extra = ["--hic"];
    assert_eq!(
        strip_sam_header(&map_stdout(c.to_str().unwrap(), &extra, &dir.join("c"), &[&r1, &r2])),
        strip_sam_header(&map_stdout(rust_bin(), &extra, &dir.join("r"), &[&r1, &r2])),
        "map --hic SAM 不一致"
    );
}

/// `--adap no` 关闭按读长插值；短读上与 C 版一致。
#[test]
fn map_adap_no_parity() {
    let Some(c) = c_minibwa() else { return };
    let Some(fa) = require_chrm("chrM-human.fa.gz") else {
        return;
    };
    let Some(reads) = require_chrm("chrM-read_1.fa.gz") else {
        return;
    };
    let dir = std::env::temp_dir().join("minibwa-rust-parity-adap-no");
    index_chrm(&c, &dir, &fa);
    let extra = ["--adap", "no"];
    assert_eq!(
        strip_sam_header(&map_stdout(c.to_str().unwrap(), &extra, &dir.join("c"), &[&reads])),
        strip_sam_header(&map_stdout(rust_bin(), &extra, &dir.join("r"), &[&reads])),
        "map --adap no SAM 不一致"
    );
}
