//! cs:Z / ds:Z 标签，对照 C 版 `cs.c` 的 `mb_write_cs_ds`。

use crate::hits::Hit;
use crate::ksw2::{CIGAR_DEL, CIGAR_INS, CIGAR_MATCH};

const CIGAR_N_SKIP: u32 = 3;
const CIGAR_EQ: u32 = 7;
const CIGAR_X: u32 = 8;
const NT: &[u8] = b"acgtn";

fn nt(b: u8) -> char {
    NT[b.min(4) as usize] as char
}

fn is_match_op(op: u32) -> bool {
    op == CIGAR_MATCH || op == CIGAR_EQ || op == CIGAR_X
}

fn write_indel_ds(out: &mut String, len: i32, seq: &[u8], ll: i32, lr: i32) {
    if ll + lr >= len {
        out.push('[');
        for &b in seq.iter().take(len as usize) {
            out.push(nt(b));
        }
        out.push(']');
        return;
    }
    let mut k = 0usize;
    if ll > 0 {
        out.push('[');
        for &b in seq[k..k + ll as usize].iter() {
            out.push(nt(b));
        }
        out.push(']');
        k += ll as usize;
    }
    let mid = (len - lr - ll) as usize;
    for &b in seq[k..k + mid].iter() {
        out.push(nt(b));
    }
    k += mid;
    if lr > 0 {
        out.push('[');
        for &b in seq[k..k + lr as usize].iter() {
            out.push(nt(b));
        }
        out.push(']');
    }
}

/// 生成 `cs:Z:` 或 `ds:Z:`（对应 `mb_write_cs_ds`）。qseq/tseq 为 nt4，与 CIGAR 对齐窗口一致。
pub fn write_cs_ds(tseq: &[u8], qseq: &[u8], r: &Hit, is_ds: bool) -> String {
    let p = r.p.as_ref().unwrap();
    let mut s = String::from(if is_ds { "ds:Z:" } else { "cs:Z:" });
    let mut q_len = 0i32;
    let mut t_len = 0i32;
    for &c in &p.cigar {
        let op = c & 0xf;
        let len = (c >> 4) as i32;
        if is_match_op(op) {
            q_len += len;
            t_len += len;
        } else if op == CIGAR_INS {
            q_len += len;
        } else if op == CIGAR_DEL || op == CIGAR_N_SKIP {
            t_len += len;
        }
    }
    let mut q_off = 0i32;
    let mut t_off = 0i32;
    for &c in &p.cigar {
        let op = c & 0xf;
        let len = (c >> 4) as i32;
        if is_match_op(op) {
            let mut l_tmp = 0i32;
            for j in 0..len {
                let cq = qseq[(q_off + j) as usize];
                let ct = tseq[(t_off + j) as usize];
                if cq != ct {
                    if l_tmp > 0 {
                        s.push_str(&format!(":{l_tmp}"));
                        l_tmp = 0;
                    }
                    s.push_str(&format!("*{}{}", nt(ct), nt(cq)));
                } else {
                    l_tmp += 1;
                }
            }
            if l_tmp > 0 {
                s.push_str(&format!(":{l_tmp}"));
            }
            q_off += len;
            t_off += len;
        } else if op == CIGAR_INS {
            if is_ds {
                let y = q_off;
                let mut z = 1i32;
                while z <= len {
                    if y - z < 0 || qseq[(y + len - z) as usize] != qseq[(y - z) as usize] {
                        break;
                    }
                    z += 1;
                }
                let lr = z - 1;
                z = 0;
                while z < len {
                    if y + len + z >= q_len || qseq[(y + len + z) as usize] != qseq[(y + z) as usize] {
                        break;
                    }
                    z += 1;
                }
                s.push('+');
                write_indel_ds(&mut s, len, &qseq[y as usize..], z, lr);
            } else {
                s.push('+');
                for j in 0..len {
                    s.push(nt(qseq[(q_off + j) as usize]));
                }
            }
            q_off += len;
        } else if op == CIGAR_DEL {
            if is_ds {
                let x = t_off;
                let mut z = 1i32;
                while z <= len {
                    if x - z < 0 || tseq[(x + len - z) as usize] != tseq[(x - z) as usize] {
                        break;
                    }
                    z += 1;
                }
                let lr = z - 1;
                z = 0;
                while z < len {
                    if x + len + z >= t_len || tseq[(x + z) as usize] != tseq[(x + len + z) as usize] {
                        break;
                    }
                    z += 1;
                }
                s.push('-');
                write_indel_ds(&mut s, len, &tseq[x as usize..], z, lr);
            } else {
                s.push('-');
                for j in 0..len {
                    s.push(nt(tseq[(t_off + j) as usize]));
                }
            }
            t_off += len;
        } else {
            // intron N
            s.push_str(&format!(
                "~{}{}{}{}{}",
                nt(tseq[t_off as usize]),
                nt(tseq[(t_off + 1) as usize]),
                len,
                nt(tseq[(t_off + len - 2) as usize]),
                nt(tseq[(t_off + len - 1) as usize])
            ));
            t_off += len;
        }
    }
    s
}
