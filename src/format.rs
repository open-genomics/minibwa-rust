//! 输出格式化，逐行对照 C 版 `format.c`。

use crate::hits::Hit;
use crate::l2bit::L2b;
use crate::opt::Opt;

/// PAF 一行（对应 `mb_fmt_paf`）。hit=None 输出未比对行。
pub fn fmt_paf(out: &mut String, l2b: &L2b, qname: &str, qlen: usize, hit: Option<&Hit>, n_seg: usize, seg_idx: usize) {
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
    if let Some(ep) = &p.p {
        let nm = p.blen - p.mlen + ep.n_ambi as i32;
        out.push_str(&format!(
            "\tNM:i:{}\tAS:i:{}\tms:i:{}\tmd:i:{}",
            nm,
            ep.dp_score,
            ep.dp_max0,
            ep.dp_max - ep.dp_max2
        ));
        if !ep.cigar.is_empty() {
            out.push_str("\tcg:Z:");
            for &c in &ep.cigar {
                out.push_str(&format!("{}{}", c >> 4, SAM_CIGAR_STR[(c & 0xf) as usize] as char));
            }
        }
        if let Some(cs) = &ep.cs {
            out.push_str(&format!("\t{}", cs));
        }
    }
    out.push('\n');
}

/// 转义 read group 行中的 `\t`/`\n`/`\r`/`\\`（对应 `mb_escape`）。
fn escape(s: &str) -> String {
    let mut out = String::new();
    let mut it = s.chars();
    while let Some(c) = it.next() {
        if c == '\\' {
            match it.next() {
                Some('t') => out.push('\t'),
                Some('n') => out.push('\n'),
                Some('r') => out.push('\r'),
                Some('\\') => out.push('\\'),
                Some(o) => {
                    out.push('\\');
                    out.push(o);
                }
                None => break,
            }
        } else {
            out.push(c);
        }
    }
    out
}

/// SAM header（对应 `mb_fmt_sam_hdr`）。返回 (header, RG 的 ID)。
pub fn fmt_sam_hdr(l2b: &L2b, ver: &str, rg: Option<&str>) -> (String, Option<String>) {
    let mut s = String::new();
    s.push_str("@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
    for ctg in &l2b.ctg {
        s.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", ctg.name, ctg.len));
    }
    let mut rg_id = None;
    if let Some(rg) = rg {
        if rg.starts_with("@RG") && !rg.contains('\t') {
            let escaped = escape(rg);
            if let Some(p) = escaped.find("\tID:") {
                let id_st = p + 4;
                let id_en = escaped[id_st..].find(['\t', '\n']).map(|i| id_st + i).unwrap_or(escaped.len());
                rg_id = Some(escaped[id_st..id_en].to_string());
                s.push_str(&escaped);
                s.push('\n');
            }
        }
    }
    s.push_str(&format!("@PG\tID:minibwa-rust\tPN:minibwa-rust\tVN:{}\n", ver));
    (s, rg_id)
}

/// CIGAR 操作符（对应 C 版 `MB_CIGAR_STR`）。
const SAM_CIGAR_STR: &[u8; 10] = b"MIDNSHP=XB";

/// 写 SEQ（反向互补如需；对应 C 版 `sam_write_sq` 的 comp=1）。
fn sam_write_sq(out: &mut String, seq: &[u8], rev: bool) {
    if rev {
        for &b in seq.iter().rev() {
            out.push(crate::dna::COMP[b as usize] as char);
        }
    } else {
        out.push_str(std::str::from_utf8(seq).unwrap_or(""));
    }
}

/// 写 CIGAR（对应 `write_sam_cigar` 的非 in_tag 分支）。
fn write_cigar(out: &mut String, qlen: i32, r: &Hit, clip_char: char) {
    if let Some(p) = &r.p {
        let clip0 = if r.rev { qlen - r.qe } else { r.qs };
        let clip1 = if r.rev { r.qs } else { qlen - r.qe };
        if clip0 > 0 {
            out.push_str(&format!("{}{}", clip0, clip_char));
        }
        for &c in &p.cigar {
            out.push_str(&format!("{}{}", c >> 4, SAM_CIGAR_STR[(c & 0xf) as usize] as char));
        }
        if clip1 > 0 {
            out.push_str(&format!("{}{}", clip1, clip_char));
        }
    } else {
        out.push('*');
    }
}

/// 找 primary（对应 `get_sam_pri`）。
pub fn get_sam_pri(hits: &[Hit]) -> Option<&Hit> {
    hits.iter().find(|h| h.sam_pri)
}

/// SAM 一行（对应 `mb_fmt_sam`，支持 SE/PE）。
///
/// `n_seg`/`seg_idx` 描述 fragment 结构；`mate_pri` 为 mate 的 primary（PE 时传入）；
/// `mate_qlen` 用于 MC:Z；`rg_id` 用于 RG:Z 标签。
#[allow(clippy::too_many_arguments)]
pub fn fmt_sam(
    out: &mut String,
    l2b: &L2b,
    qname: &str,
    seq_raw: &[u8], // ASCII 序列
    qual: Option<&[u8]>,
    hits: &[Hit],
    hit_idx: Option<usize>,
    n_seg: usize,
    seg_idx: usize,
    mate_pri: Option<&Hit>,
    mate_qlen: usize,
    rg_id: Option<&str>,
    opt: &Opt,
) {
    let r = hit_idx.and_then(|i| hits.get(i));
    let mate = mate_pri;
    // QNAME
    out.push_str(qname);
    // FLAG
    let mut flag = 0u32;
    if n_seg > 1 {
        flag |= 0x1;
    }
    if let Some(r) = r {
        if r.rev {
            flag |= 0x10;
        }
        if r.parent != r.id {
            flag |= 0x100;
        } else if !r.sam_pri {
            flag |= 0x800;
        }
    } else {
        flag |= 0x4;
    }
    if n_seg > 1 {
        if let Some(r) = r {
            if r.proper_pair {
                flag |= 0x2;
            }
        }
        if seg_idx == 0 {
            flag |= 0x40;
        } else if seg_idx == n_seg - 1 {
            flag |= 0x80;
        }
        if mate.is_none() {
            flag |= 0x8;
        } else if let Some(m) = mate {
            if m.rev {
                flag |= 0x20;
            }
        }
    }
    out.push_str(&format!("\t{}", flag));

    // RNAME POS MAPQ CIGAR
    let mut this_tid: i64 = -1;
    let mut this_pos: i64 = -1;
    match r {
        None => {
            if let Some(m) = mate {
                this_tid = m.tid;
                this_pos = m.ts;
                out.push_str(&format!("\t{}\t{}\t0\t*", l2b.ctg[m.tid as usize].name, m.ts + 1));
            } else {
                out.push_str("\t*\t0\t0\t*");
            }
        }
        Some(r) => {
            this_tid = r.tid;
            this_pos = r.ts;
            out.push_str(&format!("\t{}\t{}\t{}", l2b.ctg[r.tid as usize].name, r.ts + 1, r.mapq));
            out.push('\t');
            if r.p.is_some() {
                let clip_char = if flag & 0x800 != 0 { 'H' } else { 'S' };
                write_cigar(out, seq_raw.len() as i32, r, clip_char);
            } else {
                out.push('*');
            }
        }
    }

    // RNEXT PNEXT TLEN
    if n_seg > 1 {
        let mut tlen = 0i64;
        if let Some(m) = mate {
            if this_tid >= 0 {
                if this_tid == m.tid {
                    if let Some(rr) = r {
                        let this_pos5 = if rr.rev { rr.te - 1 } else { this_pos };
                        let next_pos5 = if m.rev { m.te - 1 } else { m.ts };
                        tlen = next_pos5 - this_pos5;
                    }
                    out.push_str("\t=\t");
                } else {
                    out.push_str(&format!("\t{}\t", l2b.ctg[m.tid as usize].name));
                }
                out.push_str(&format!("{}\t", m.ts + 1));
            } else {
                // this_tid < 0：本段未比对，放在 mate 位置
                out.push_str(&format!("\t{}\t{}\t", l2b.ctg[m.tid as usize].name, m.ts + 1));
            }
        } else if this_tid >= 0 {
            // mate 未比对：放在本段 primary 处
            let r_pri = get_sam_pri(hits).expect("有 hit 必有 sam_pri");
            let cur = r.expect("this_tid>=0 必有 hit");
            if cur.tid == r_pri.tid {
                out.push_str(&format!("\t=\t{}\t", r_pri.ts + 1));
            } else {
                out.push_str(&format!("\t{}\t{}\t", l2b.ctg[r_pri.tid as usize].name, r_pri.ts + 1));
            }
        } else {
            out.push_str("\t*\t0\t");
        }
        if tlen > 0 {
            tlen += 1;
        } else if tlen < 0 {
            tlen -= 1;
        }
        out.push_str(&format!("{}\t", tlen));
    } else {
        out.push_str("\t*\t0\t0\t");
    }

    // SEQ QUAL
    match r {
        None => {
            sam_write_sq(out, seq_raw, false);
            out.push('\t');
            match qual {
                Some(q) => out.push_str(std::str::from_utf8(q).unwrap_or("*")),
                None => out.push('*'),
            }
        }
        Some(r) => {
            if flag & 0x900 == 0 {
                sam_write_sq(out, seq_raw, r.rev);
            } else if flag & 0x100 != 0 {
                out.push('*');
            } else {
                // 补充比对：片段 SEQ
                let seg = &seq_raw[r.qs as usize..r.qe as usize];
                sam_write_sq(out, seg, r.rev);
            }
            out.push('\t');
            match qual {
                Some(q) if flag & 0x900 == 0 => {
                    let qseg = if r.rev {
                        q.iter().rev().copied().collect::<Vec<_>>()
                    } else {
                        q.to_vec()
                    };
                    out.push_str(std::str::from_utf8(&qseg).unwrap_or("*"));
                }
                Some(_) if flag & 0x100 != 0 => out.push('*'),
                Some(q) => {
                    // 补充比对：片段 QUAL（反转不互补）
                    let qseg = &q[r.qs as usize..r.qe as usize];
                    if r.rev {
                        for &b in qseg.iter().rev() {
                            out.push(b as char);
                        }
                    } else {
                        out.push_str(std::str::from_utf8(qseg).unwrap_or("*"));
                    }
                }
                None => out.push('*'),
            }
        }
    }

    // tags
    if let Some(r) = r {
        if let Some(rg) = rg_id {
            out.push_str(&format!("\tRG:Z:{}", rg));
        }
        if let Some(p) = &r.p {
            let nm = r.blen - r.mlen + p.n_ambi as i32;
            out.push_str(&format!(
                "\tNM:i:{}\tAS:i:{}\tms:i:{}\tmd:i:{}",
                nm,
                p.dp_score,
                p.dp_max0,
                p.dp_max - p.dp_max2
            ));
            // MC:Z mate CIGAR 与 MQ:i mate MAPQ（r_next 即 mate 的 primary）
            if n_seg > 1 {
                if let Some(m) = mate {
                    if m.p.is_some() && !m.p.as_ref().unwrap().cigar.is_empty() && mate_qlen > 0 {
                        out.push_str("\tMC:Z:");
                        write_cigar(out, mate_qlen as i32, m, 'S');
                        out.push_str(&format!("\tMQ:i:{}", m.mapq));
                    }
                }
            }
            // cs/ds/MD 标签（对应 C 版 write_tags 之后的 cs 输出）
            if let Some(cs) = &p.cs {
                out.push_str(&format!("\t{}", cs));
            }
            // SA:Z（对应 C 版：parent==id 且有其他 primary hit 时列出其余 primary）
            if r.parent == r.id && hits.len() > 1 {
                let hi = hit_idx.unwrap();
                let n_sa = hits
                    .iter()
                    .enumerate()
                    .filter(|(i, q)| *i != hi && q.parent == q.id && q.p.is_some())
                    .count();
                if n_sa > 0 {
                    out.push_str("\tSA:Z:");
                    let qlen = seq_raw.len() as i64;
                    for (i, q) in hits.iter().enumerate() {
                        if i == hi || q.parent != q.id || q.p.is_none() {
                            continue;
                        }
                        let (qs, qe) = (q.qs as i64, q.qe as i64);
                        let (ts, te) = (q.ts, q.te);
                        let (l_m, l_i, l_d) = if qe - qs < te - ts {
                            let l_m = qe - qs;
                            (l_m, 0, te - ts - l_m)
                        } else {
                            let l_m = te - ts;
                            (l_m, qe - qs - l_m, 0)
                        };
                        let clip5 = if q.rev { qlen - qe } else { qs };
                        let clip3 = if q.rev { qs } else { qlen - qe };
                        out.push_str(&format!(
                            "{},{},{},",
                            l2b.ctg[q.tid as usize].name,
                            ts + 1,
                            if q.rev { '-' } else { '+' }
                        ));
                        if clip5 > 0 {
                            out.push_str(&format!("{}S", clip5));
                        }
                        if l_m > 0 {
                            out.push_str(&format!("{}M", l_m));
                        }
                        if l_i > 0 {
                            out.push_str(&format!("{}I", l_i));
                        }
                        if l_d > 0 {
                            out.push_str(&format!("{}D", l_d));
                        }
                        if clip3 > 0 {
                            out.push_str(&format!("{}S", clip3));
                        }
                        out.push_str(&format!(
                            ",{},{};",
                            q.mapq,
                            q.blen - q.mlen + q.p.as_ref().unwrap().n_ambi as i32
                        ));
                    }
                }
                // n2:i / XA:Z（对应 C 版 write_tags 后的 xa_max 分支）
                if opt.xa_max > 0 {
                    let dp_best = r.p.as_ref().map_or(0, |pp| pp.dp_max);
                    let mut n_xa = 0;
                    for (i, q) in hits.iter().enumerate() {
                        let Some(qp) = q.p.as_ref() else { continue };
                        if i != hi && q.parent as usize == hi && qp.dp_max as f64 >= opt.out_s as f64 * dp_best as f64 {
                            n_xa += 1;
                        }
                    }
                    if n_xa > 0 {
                        out.push_str(&format!("\tn2:i:{}", n_xa));
                    }
                    if n_xa > 0 && n_xa <= opt.xa_max {
                        out.push_str("\tXA:Z:");
                        for (i, q) in hits.iter().enumerate() {
                            let Some(qp) = q.p.as_ref() else { continue };
                            if i != hi
                                && q.parent as usize == hi
                                && qp.dp_max as f64 >= opt.out_s as f64 * dp_best as f64
                            {
                                out.push_str(&format!(
                                    "{},{}{},",
                                    l2b.ctg[q.tid as usize].name,
                                    if q.rev { '-' } else { '+' },
                                    q.ts + 1
                                ));
                                write_cigar(out, seq_raw.len() as i32, q, 'S');
                                out.push_str(&format!(",{};", q.blen - q.mlen + qp.n_ambi as i32));
                            }
                        }
                    }
                }
            }
        }
    }
    out.push('\n');
}
