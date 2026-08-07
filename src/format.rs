//! 输出格式化，逐行对照 C 版 `format.c`。

use crate::hits::Hit;
use crate::l2bit::L2b;

/// PAF 一行（对应 `mb_fmt_paf`）。hit=None 输出未比对行。
pub fn fmt_paf(
    out: &mut String,
    l2b: &L2b,
    qname: &str,
    qlen: usize,
    hit: Option<&Hit>,
    n_seg: usize,
    seg_idx: usize,
) {
    out.push_str(qname);
    if n_seg > 1 {
        out.push_str(&format!("/{}", seg_idx + 1));
    }
    out.push_str(&format!("\t{qlen}"));
    let Some(p) = hit else {
        out.push_str("\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n");
        return;
    };
    let ctg = &l2b.ctg[p.tid as usize];
    out.push_str(&format!(
        "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttp:A:{}\ts1:i:{}\tcm:i:{}",
        p.qs,
        p.qe,
        if p.rev { '-' } else { '+' },
        ctg.name,
        ctg.len,
        p.ts,
        p.te,
        p.mlen,
        p.blen,
        p.mapq,
        if p.parent == p.id { 'P' } else { 'S' },
        p.score,
        p.cnt
    ));
    if p.parent == p.id {
        out.push_str(&format!("\ts2:i:{}", p.subsc.max(0)));
    }
    // M5：dp 标签与 cg/cs 将在此补充
    out.push('\n');
}

/// SAM header（对应 `mb_fmt_sam_hdr`）。
pub fn fmt_sam_hdr(l2b: &L2b, ver: &str) -> String {
    let mut s = String::new();
    s.push_str("@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
    for ctg in &l2b.ctg {
        s.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", ctg.name, ctg.len));
    }
    s.push_str(&format!("@PG\tID:minibwa-rust\tPN:minibwa-rust\tVN:{}\n", ver));
    s
}

const SAM_CIGAR_STR: &[u8; 9] = b"MIDNSHP=X";

/// 写 SEQ（反向互补如需）。
fn sam_write_sq(out: &mut String, seq: &[u8], rev: bool) {
    if rev {
        for &b in seq.iter().rev() {
            out.push(crate::dna::COMP[b as usize] as char);
        }
    } else {
        out.push_str(std::str::from_utf8(seq).unwrap_or(""));
    }
}

/// SAM 一行（对应 `mb_fmt_sam`，SE 简化版）。
pub fn fmt_sam(
    out: &mut String,
    l2b: &L2b,
    qname: &str,
    seq_raw: &[u8], // ASCII 序列
    qual: Option<&[u8]>,
    hit: Option<&Hit>,
) {
    // QNAME
    out.push_str(qname);
    // FLAG
    let mut flag = 0u32;
    let r = hit;
    if r.is_none() {
        flag |= 0x4;
    } else {
        let r = r.unwrap();
        if r.rev { flag |= 0x10; }
        if r.parent != r.id { flag |= 0x100; }
        else if !r.sam_pri { flag |= 0x800; }
    }
    out.push_str(&format!("\t{}", flag));

    // RNAME POS MAPQ CIGAR
    match r {
        None => out.push_str("\t*\t0\t0\t*"),
        Some(r) => {
            out.push_str(&format!("\t{}\t{}\t{}", l2b.ctg[r.tid as usize].name, r.ts + 1, r.mapq));
            // CIGAR
            if let Some(p) = &r.p {
                let clip0 = if r.rev { seq_raw.len() as i32 - r.qe } else { r.qs };
                let clip1 = if r.rev { r.qs } else { seq_raw.len() as i32 - r.qe };
                let clip_char = if flag & 0x900 != 0 && !(false) { 'H' } else { 'S' }; // SUPP_SOFT 未支持
                out.push('\t');
                if clip0 > 0 { out.push_str(&format!("{}{}", clip0, clip_char)); }
                for &c in &p.cigar {
                    out.push_str(&format!("{}{}", c >> 4, SAM_CIGAR_STR[(c & 0xf) as usize] as char));
                }
                if clip1 > 0 { out.push_str(&format!("{}{}", clip1, clip_char)); }
            } else {
                out.push_str("\t*");
            }
        }
    }

    // RNEXT PNEXT TLEN（SE）
    out.push_str("\t*\t0\t0\t");

    // SEQ QUAL
    match r {
        None => {
            sam_write_sq(out, seq_raw, false);
            out.push('\t');
            match qual { Some(q) => out.push_str(std::str::from_utf8(q).unwrap_or("*")), None => out.push('*') }
        }
        Some(r) => {
            if flag & 0x900 == 0 {
                sam_write_sq(out, seq_raw, r.rev);
            } else if flag & 0x100 != 0 {
                // 次要比对：可能省略 SEQ
                out.push('*');
            } else {
                // 补充比对
                let seg = &seq_raw[r.qs as usize..r.qe as usize];
                sam_write_sq(out, seg, r.rev);
            }
            out.push('\t');
            match qual {
                Some(q) if flag & 0x900 == 0 => { let qseg = if r.rev { q.iter().rev().copied().collect::<Vec<_>>() } else { q.to_vec() }; out.push_str(std::str::from_utf8(&qseg).unwrap_or("*")); }
                _ => out.push('*'),
            }
        }
    }

    // tags
    if let Some(r) = r {
        if let Some(p) = &r.p {
            let nm = r.blen - r.mlen + p.n_ambi as i32;
            out.push_str(&format!("\tNM:i:{}\tAS:i:{}\tms:i:{}\tmd:i:{}", nm, p.dp_score, p.dp_max0, p.dp_max - p.dp_max2));
        }
    }
    out.push('\n');
}
