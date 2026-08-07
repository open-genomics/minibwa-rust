use anyhow::Result;
use std::io::Write;

/// SAM flag constants
pub mod flags {
    /// Read paired
    pub const PAIRED: u16 = 0x1;
    /// Read mapped in proper pair
    pub const PROPER_PAIR: u16 = 0x2;
    /// Read unmapped
    pub const UNMAP: u16 = 0x4;
    /// Mate unmapped
    pub const MUNMAP: u16 = 0x8;
    /// Read reverse strand
    pub const REVERSE: u16 = 0x10;
    /// Mate reverse strand
    pub const MREVERSE: u16 = 0x20;
    /// First in pair
    pub const READ1: u16 = 0x40;
    /// Second in pair
    pub const READ2: u16 = 0x80;
    /// Secondary alignment
    pub const SECONDARY: u16 = 0x100;
    /// QC failure
    pub const QCFAIL: u16 = 0x200;
    /// Duplicate
    pub const DUP: u16 = 0x400;
    /// Supplementary alignment
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// Write SAM header (@HD, @SQ, @PG) to output
pub fn write_header<W: Write, S: AsRef<str>>(out: &mut W, contigs: &[(S, u32)]) -> Result<()> {
    writeln!(out, "@HD\tVN:1.6\tSO:unsorted")?;
    for (name, len) in contigs {
        writeln!(out, "@SQ\tSN:{}\tLN:{}", name.as_ref(), len)?;
    }
    writeln!(
        out,
        "@PG\tID:minibwa-rust\tPN:minibwa-rust\tVN:{}",
        env!("CARGO_PKG_VERSION")
    )?;
    Ok(())
}

/// Format an unmapped SAM record (FLAG=4)
pub fn format_unmapped(qname: &str, seq: &str, qual: &str) -> String {
    format!("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}", qname, seq, qual,)
}

/// Format a mapped SAM record with optional tags
pub fn format_record(
    qname: &str,
    flag: u16,
    rname: &str,
    pos: u32,
    mapq: u8,
    cigar: &str,
    seq: &str,
    qual: &str,
    score: i32,
    sub_score: i32,
    nm: u32,
) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tAS:i:{}\tXS:i:{}\tNM:i:{}",
        qname, flag, rname, pos, mapq, cigar, seq, qual, score, sub_score, nm,
    )
}

/// Format a mapped SAM record with MD:Z and SA:Z tags
pub fn format_record_with_md_sa(
    qname: &str,
    flag: u16,
    rname: &str,
    pos: u32,
    mapq: u8,
    cigar: &str,
    seq: &str,
    qual: &str,
    score: i32,
    sub_score: i32,
    nm: u32,
    md_tag: &str,
    sa_tag: &str,
) -> String {
    if sa_tag.is_empty() {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tAS:i:{}\tXS:i:{}\tNM:i:{}\tMD:Z:{}",
            qname, flag, rname, pos, mapq, cigar, seq, qual, score, sub_score, nm, md_tag,
        )
    } else {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tAS:i:{}\tXS:i:{}\tNM:i:{}\tMD:Z:{}\tSA:Z:{}",
            qname, flag, rname, pos, mapq, cigar, seq, qual, score, sub_score, nm, md_tag, sa_tag,
        )
    }
}

pub fn format_record_with_optional_tags(
    qname: &str,
    flag: u16,
    rname: &str,
    pos: u32,
    mapq: u8,
    cigar: &str,
    seq: &str,
    qual: &str,
    score: i32,
    sub_score: i32,
    nm: u32,
    md_tag: &str,
    sa_tag: &str,
) -> String {
    let mut line = format_record(qname, flag, rname, pos, mapq, cigar, seq, qual, score, sub_score, nm);
    if !md_tag.is_empty() {
        line.push_str("\tMD:Z:");
        line.push_str(md_tag);
    }
    if !sa_tag.is_empty() {
        line.push_str("\tSA:Z:");
        line.push_str(sa_tag);
    }
    line
}

/// Generate MD:Z tag from reference and query sequences aligned according to CIGAR.
///
/// The MD:Z tag encodes the reference sequence at mismatch positions for variant calling.
/// Format: numbers for matches, letters for mismatches, ^ followed by letters for deletions.
///
/// # Arguments
/// * `reference` - Reference sequence segment (already oriented to match query strand)
/// * `query` - Query sequence segment (oriented to match the alignment strand)
/// * `cigar` - CIGAR string describing the alignment
///
/// # Returns
/// MD:Z tag string
///
/// # Example
/// ```
/// // Perfect match: all bases align
/// let md = minibwa_rust::io::sam::generate_md_tag(b"ACGTACGT", b"ACGTACGT", "8M");
/// assert_eq!(md, "8");
///
/// // One mismatch: ref A at pos 4, query T
/// let md = minibwa_rust::io::sam::generate_md_tag(b"ACGTACGT", b"ACGTTCGT", "8M");
/// assert_eq!(md, "4A3");
/// ```
pub fn generate_md_tag(reference: &[u8], query: &[u8], cigar: &str) -> String {
    let ops = parse_cigar_ops(cigar);
    let mut md = String::new();
    let mut ref_pos = 0usize;
    let mut query_pos = 0usize;
    let mut match_count = 0usize;

    for (op, len) in ops {
        match op {
            'M' | '=' | 'X' => {
                // Alignment match: compare query vs reference
                for _ in 0..len {
                    if ref_pos >= reference.len() || query_pos >= query.len() {
                        break;
                    }
                    let ref_base = reference[ref_pos].to_ascii_uppercase();
                    let query_base = query[query_pos].to_ascii_uppercase();

                    if ref_base == query_base {
                        match_count += 1;
                    } else {
                        // Mismatch: output accumulated matches then the mismatched ref base
                        if match_count > 0 {
                            md.push_str(&match_count.to_string());
                            match_count = 0;
                        }
                        md.push(ref_base as char);
                    }
                    ref_pos += 1;
                    query_pos += 1;
                }
            }
            'I' => {
                // Insertion: skip query bases, not in MD tag
                query_pos += len;
            }
            'D' | 'N' => {
                // Deletion/skip: output ^ followed by deleted reference bases
                if match_count > 0 {
                    md.push_str(&match_count.to_string());
                    match_count = 0;
                }
                md.push('^');
                for _ in 0..len {
                    if ref_pos < reference.len() {
                        md.push(reference[ref_pos].to_ascii_uppercase() as char);
                        ref_pos += 1;
                    }
                }
            }
            'S' => {
                // Soft clip: skip query bases, not in MD tag
                query_pos += len;
            }
            'H' => {
                // Hard clip: nothing consumed from either sequence
            }
            'P' => {
                // Padding: skip both (rare)
                ref_pos += len;
                query_pos += len;
            }
            _ => {
                // Unknown operator: skip
            }
        }
    }

    // Output any remaining match count
    if match_count > 0 {
        md.push_str(&match_count.to_string());
    }

    md
}

/// Parse CIGAR string into (operator, length) pairs.
fn parse_cigar_ops(cigar: &str) -> Vec<(char, usize)> {
    let mut result = Vec::new();
    let mut num = 0usize;

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch as usize - '0' as usize);
        } else if num > 0 {
            result.push((ch, num));
            num = 0;
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_format() {
        let mut buf = Vec::new();
        let contigs = vec![("chr1", 1000u32), ("chr2", 2000u32)];
        write_header(&mut buf, &contigs).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("@HD\tVN:1.6\tSO:unsorted"));
        assert!(s.contains("@SQ\tSN:chr1\tLN:1000"));
        assert!(s.contains("@SQ\tSN:chr2\tLN:2000"));
        assert!(s.contains("@PG\tID:minibwa-rust"));
    }

    #[test]
    fn unmapped_format() {
        let line = format_unmapped("read1", "ACGT", "IIII");
        assert!(line.contains("\t4\t"));
        assert!(line.starts_with("read1\t"));
    }

    #[test]
    fn record_format() {
        let line = format_record("read1", 0, "chr1", 100, 60, "50M", "ACGT", "IIII", 100, 0, 2);
        assert!(line.starts_with("read1\t0\tchr1\t100\t60\t50M\t"));
        assert!(line.contains("AS:i:100"));
        assert!(line.contains("NM:i:2"));
    }

    #[test]
    fn header_empty_contigs() {
        let mut buf = Vec::new();
        let contigs: Vec<(&str, u32)> = vec![];
        write_header(&mut buf, &contigs).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("@HD"));
        assert!(s.contains("@PG"));
        assert!(!s.contains("@SQ"));
    }

    #[test]
    fn unmapped_has_correct_tab_count() {
        let line = format_unmapped("r1", "ACGT", "IIII");
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 11);
        assert_eq!(fields[0], "r1");
        assert_eq!(fields[1], "4");
        assert_eq!(fields[2], "*");
        assert_eq!(fields[3], "0");
        assert_eq!(fields[9], "ACGT");
        assert_eq!(fields[10], "IIII");
    }

    #[test]
    fn record_format_reverse_complement() {
        let line = format_record("read1", 16, "chr1", 50, 30, "20M", "ACGT", "IIII", 40, 10, 1);
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[1], "16");
        assert!(line.contains("XS:i:10"));
    }

    #[test]
    fn record_format_secondary_alignment() {
        let line = format_record("read1", 256, "chr2", 200, 0, "10M1I10M", "ACGT", "IIII", 30, 50, 3);
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[1], "256");
        assert_eq!(fields[2], "chr2");
        assert_eq!(fields[3], "200");
        assert_eq!(fields[4], "0");
        assert_eq!(fields[5], "10M1I10M");
    }

    #[test]
    fn header_with_string_contigs() {
        let mut buf = Vec::new();
        let contigs = vec![("chrX".to_string(), 155_270_560u32)];
        write_header(&mut buf, &contigs).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("@SQ\tSN:chrX\tLN:155270560"));
    }

    #[test]
    fn md_tag_exact_match() {
        // Perfect match: all bases match
        let ref_seq = b"ACGTACGTACGT";
        let query = b"ACGTACGTACGT";
        let md = generate_md_tag(ref_seq, query, "12M");
        assert_eq!(md, "12");
    }

    #[test]
    fn md_tag_single_mismatch() {
        // One mismatch in the middle
        let ref_seq = b"ACGTACGTACGT";
        let query = b"ACGTTCGTACGT"; // T at position 4 instead of A
        let md = generate_md_tag(ref_seq, query, "12M");
        assert_eq!(md, "4A7");
    }

    #[test]
    fn md_tag_multiple_mismatches() {
        // Multiple mismatches
        let ref_seq = b"ACGTACGTACGT";
        let query = b"TCGTTCGTACGT"; // T at pos 0, T at pos 4
        let md = generate_md_tag(ref_seq, query, "12M");
        assert_eq!(md, "A3A7");
    }

    #[test]
    fn md_tag_with_deletion() {
        // Deletion in reference: CIGAR 4M3D4M
        // Query consumes: 4 + 0 + 4 = 8 bases
        // Ref consumes: 4 + 3 + 4 = 11 bases
        let ref_seq = b"ACGTTTTACGT"; // 11 bases: ACGT + TTT + ACGT
        let query = b"ACGTACGT"; // 8 bases
        let md = generate_md_tag(ref_seq, query, "4M3D4M");
        assert_eq!(md, "4^TTT4");
    }

    #[test]
    fn md_tag_with_insertion() {
        // Insertion in query (not in MD tag)
        let ref_seq = b"ACGTACGT";
        let query = b"ACGTTACGT"; // extra T at position 4
        let md = generate_md_tag(ref_seq, query, "4M1I4M");
        assert_eq!(md, "8");
    }

    #[test]
    fn md_tag_with_soft_clip() {
        // Soft clipping at both ends
        let ref_seq = b"ACGTACGT";
        let query = b"NNACGTACGTNN";
        let md = generate_md_tag(ref_seq, query, "2S8M2S");
        assert_eq!(md, "8");
    }

    #[test]
    fn md_tag_complex() {
        // Simple mismatch test
        // CIGAR: 4M1X4M
        // Ref consumes: 4 + 1 + 4 = 9 bases
        // Query consumes: 4 + 1 + 4 = 9 bases
        let ref_seq = b"ACGTAACGT"; // 9 bases
        let query = b"ACGTTACGT"; // 9 bases (T mismatch at pos 4)
        let md = generate_md_tag(ref_seq, query, "4M1X4M");
        // 4M = ACGT matches
        // 1X = A->T mismatch, output "A"
        // 4M = ACGT matches
        assert_eq!(md, "4A4");
    }

    #[test]
    fn md_tag_with_deletion_and_mismatch() {
        // Deletion and mismatch combined
        // CIGAR: 4M3D4M
        // Query consumes: 4 + 0 + 4 = 8 bases
        // Ref consumes: 4 + 3 + 4 = 11 bases
        let ref_seq = b"ACGTTTTACGT"; // 11 bases: ACGT + TTT + ACGT
        let query = b"ACGTTCGT"; // 8 bases with T mismatch at pos 5 of alignment
                                 // 4M = ACGT matches, 3D = TTT deleted, 4M = query[4:8]=TCGT vs ref[7:11]=ACGT
                                 // First M of last 4M: ref=A, query=T -> mismatch A
                                 // Next 3 M: CGT matches
        let md = generate_md_tag(ref_seq, query, "4M3D4M");
        assert_eq!(md, "4^TTTA3");
    }

    #[test]
    fn format_record_with_md_tag() {
        let line = format_record_with_md_sa("read1", 0, "chr1", 100, 60, "50M", "ACGT", "IIII", 100, 0, 2, "50", "");
        assert!(line.contains("MD:Z:50"));
        assert!(line.contains("AS:i:100"));
    }

    #[test]
    fn format_record_with_md_sa_tag() {
        let line = format_record_with_md_sa(
            "read1",
            0,
            "chr1",
            100,
            60,
            "50M",
            "ACGT",
            "IIII",
            100,
            0,
            2,
            "50",
            "chr2,200,+,50M,60,0;",
        );
        assert!(line.contains("MD:Z:50"));
        assert!(line.contains("SA:Z:chr2,200,+,50M,60,0;"));
        assert!(line.contains("AS:i:100"));
    }

    #[test]
    fn format_record_with_sa_without_md_tag() {
        let line = format_record_with_optional_tags(
            "read1",
            0,
            "chr1",
            100,
            60,
            "50M",
            "ACGT",
            "IIII",
            100,
            0,
            2,
            "",
            "chr2,200,+,50M,60,0;",
        );
        assert!(line.contains("SA:Z:chr2,200,+,50M,60,0;"));
        assert!(!line.contains("MD:Z:"));
    }

    #[test]
    fn md_tag_case_insensitive() {
        // Mixed case should work
        let ref_seq = b"acgtACGT";
        let query = b"ACGTacgt";
        let md = generate_md_tag(ref_seq, query, "8M");
        assert_eq!(md, "8");
    }

    #[test]
    fn md_tag_empty() {
        let md = generate_md_tag(b"", b"", "");
        assert_eq!(md, "");
    }
}
