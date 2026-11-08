//! 损坏/伪造索引不得在读取前按文件头声明值分配巨量内存。
//!
//! 审计复现：`.mbw` 的 l2[4]（seq_len）改为 1<<50、`.l2b` 的 n_ctg 改为
//! 1<<60 时，旧加载器分别 abort（memory allocation failed）与 panic
//! （capacity overflow）。这里的断言只要求返回可诊断错误。

use minibwa_rust::bwt::{self, Bwt};
use minibwa_rust::l2bit::{self, L2b};
use std::fs;
use std::path::PathBuf;

fn temp_path(tag: &str) -> PathBuf {
    std::env::temp_dir().join(format!("minibwa-corrupt-{}-{}.idx", std::process::id(), tag))
}

fn write_temp(tag: &str, bytes: &[u8]) -> PathBuf {
    let path = temp_path(tag);
    fs::write(&path, bytes).expect("write temp index");
    path
}

fn push_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

/// `.mbw`: magic + sa_bit + primary + l2[1..5]（l2[4]=seq_len）+ data + n_sa + sa
fn bwt_bytes(sa_bit: u32, seq_len: u64, data_words: usize, n_sa: u64) -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend_from_slice(bwt::MAGIC);
    buf.extend_from_slice(&sa_bit.to_le_bytes());
    push_u64(&mut buf, 0); // primary
    push_u64(&mut buf, 0); // l2[1]
    push_u64(&mut buf, 0); // l2[2]
    push_u64(&mut buf, 0); // l2[3]
    push_u64(&mut buf, seq_len); // l2[4]
    for _ in 0..data_words {
        push_u64(&mut buf, 0);
    }
    push_u64(&mut buf, n_sa);
    buf
}

/// `.l2b`: magic + dummy + 7 个 u64 头部字段。
fn l2b_header(n_ctg: u64, tot_len: u64, n_ambi: u64, n_mask: u64, len_name: u64, len_comm: u64, n_pac: u64) -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend_from_slice(l2bit::MAGIC);
    buf.extend_from_slice(&0u32.to_le_bytes());
    push_u64(&mut buf, n_ctg);
    push_u64(&mut buf, tot_len);
    push_u64(&mut buf, n_ambi);
    push_u64(&mut buf, n_mask);
    push_u64(&mut buf, len_name);
    push_u64(&mut buf, len_comm);
    push_u64(&mut buf, n_pac);
    buf
}

#[test]
fn bwt_load_rejects_huge_seq_len_without_allocating() {
    let path = write_temp("bwt-huge-seq", &bwt_bytes(u32::MAX, 1 << 50, 0, 0));
    let err = match Bwt::load(&path) {
        Ok(_) => panic!("expected Bwt::load to reject the forged index"),
        Err(e) => e,
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("剩余字节"), "expected a size-mismatch error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn bwt_load_rejects_invalid_sa_bit() {
    // seq_len=64 → data_len=2 words；sa_bit=64 是非法移位值。
    let path = write_temp("bwt-bad-sa-bit", &bwt_bytes(64, 64, 12, 1));
    let err = match Bwt::load(&path) {
        Ok(_) => panic!("expected Bwt::load to reject the forged index"),
        Err(e) => e,
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("sa_bit"), "expected sa_bit error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn bwt_load_rejects_zero_sa_entries() {
    // sa_bit 有效但 n_sa=0：旧实现放行，map 阶段 cache(10) 越界。
    let path = write_temp("bwt-zero-sa", &bwt_bytes(3, 64, 12, 0));
    let err = match Bwt::load(&path) {
        Ok(_) => panic!("expected Bwt::load to reject the forged index"),
        Err(e) => e,
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("n_sa"), "expected n_sa error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn bwt_load_rejects_out_of_range_primary() {
    let mut bytes = bwt_bytes(u32::MAX, 64, 12, 0);
    bytes[8..16].copy_from_slice(&64u64.to_le_bytes()); // primary == seq_len
    let path = write_temp("bwt-bad-primary", &bytes);
    let err = match Bwt::load(&path) {
        Ok(_) => panic!("expected Bwt::load to reject the forged index"),
        Err(e) => e,
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("primary"), "expected primary error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn bwt_load_rejects_non_monotonic_l2() {
    let mut bytes = bwt_bytes(u32::MAX, 64, 12, 0);
    bytes[24..32].copy_from_slice(&100u64.to_le_bytes()); // l2[2] > l2[4]
    let path = write_temp("bwt-bad-l2", &bytes);
    let err = match Bwt::load(&path) {
        Ok(_) => panic!("expected Bwt::load to reject the forged index"),
        Err(e) => e,
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("L2"), "expected L2 error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn l2b_rejects_huge_n_ctg_without_allocating() {
    let path = write_temp("l2b-huge-ctg", &l2b_header(1 << 60, 0, 0, 0, 0, 0, 0));
    let err = L2b::from_bytes(&fs::read(&path).unwrap()).unwrap_err();
    let msg = format!("{err:#}");
    assert!(msg.contains("剩余字节"), "expected a size-mismatch error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn l2b_rejects_huge_n_pac_without_allocating() {
    let path = write_temp("l2b-huge-pac", &l2b_header(0, 1 << 40, 0, 0, 0, 0, 1 << 50));
    let err = L2b::from_bytes(&fs::read(&path).unwrap()).unwrap_err();
    let msg = format!("{err:#}");
    assert!(msg.contains("剩余字节"), "expected a size-mismatch error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn l2b_rejects_inconsistent_n_pac() {
    let path = write_temp("l2b-bad-pac", &l2b_header(0, 64, 0, 0, 0, 0, 0));
    let err = L2b::from_bytes(&fs::read(&path).unwrap()).unwrap_err();
    let msg = format!("{err:#}");
    assert!(msg.contains("n_pac"), "expected an n_pac consistency error, got: {msg}");
    let _ = fs::remove_file(path);
}

#[test]
fn l2b_load_rejects_comment_stream_shorter_than_contigs() {
    // 头部 len_comm=0 但 n_ctg=1：旧实现直接 `comms[p_comm]` 索引 panic。
    let mut buf = l2b_header(1, 4, 0, 0, 2, 0, 1);
    push_u64(&mut buf, 4); // contig len，总和 = tot_len
    push_u64(&mut buf, 0); // pac word
    buf.extend_from_slice(b"A\0");
    let path = write_temp("l2b-short-comment", &buf);

    let err = match L2b::load(&path) {
        Ok(_) => panic!("expected L2b::load to reject the forged index"),
        Err(e) => e,
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("注释串长度小于 contig 数"), "unexpected error: {msg}");
    let _ = fs::remove_file(path);
}
