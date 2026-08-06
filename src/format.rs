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
