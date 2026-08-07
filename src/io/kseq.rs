//! kseq 兼容的 FASTA/FASTQ 迭代器（对应 C 版 kseq.h 的 `kseq_read`）。
//!
//! 语义与 kseq 一致：
//! - name 取到第一个空白字符；comment 为头行剩余部分（不跳过前导空白）；
//! - 序列行跳过 '\n' 但保留其他字节；
//! - FASTQ 的 '+' 行剩余部分被跳过，qual 读到与 seq 等长。

use std::io::{BufRead, BufReader, Read};

#[derive(Debug, Clone, Default)]
pub struct SeqRecord {
    pub name: String,
    pub comment: Option<String>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

pub struct SeqIter<R: Read> {
    reader: BufReader<R>,
    pending_header: Option<Vec<u8>>,
    done: bool,
}

impl<R: Read> SeqIter<R> {
    pub fn new(reader: R) -> Self {
        SeqIter {
            reader: BufReader::new(reader),
            pending_header: None,
            done: false,
        }
    }

    fn read_line_bytes(&mut self, buf: &mut Vec<u8>) -> std::io::Result<bool> {
        buf.clear();
        let n = self.reader.read_until(b'\n', buf)?;
        if n == 0 {
            return Ok(false);
        }
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        Ok(true)
    }

    /// 读一条记录；EOF 返回 None（对应 kseq_read 返回 -1）。
    pub fn next_record(&mut self) -> Option<std::io::Result<SeqRecord>> {
        if self.done {
            return None;
        }
        // 头行内容（不含 '>'/'@' 与行尾 '\n'）
        let header: Vec<u8> = if let Some(h) = self.pending_header.take() {
            h
        } else {
            loop {
                let mut b = [0u8; 1];
                match self.reader.read(&mut b) {
                    Ok(0) => {
                        self.done = true;
                        return None;
                    }
                    Err(e) => return Some(Err(e)),
                    Ok(_) => {
                        if b[0] == b'>' || b[0] == b'@' {
                            break;
                        }
                    }
                }
            }
            // '>'/'@' 已消费，读入头行剩余部分
            let mut h = Vec::new();
            match self.read_line_bytes(&mut h) {
                Ok(false) => {
                    self.done = true;
                    return None;
                }
                Err(e) => return Some(Err(e)),
                Ok(true) => {}
            }
            h
        };

        let mut rec = SeqRecord::default();
        let ws = header.iter().position(|&b| matches!(b, b' ' | b'\t' | 0x0b | 0x0c | b'\r'));
        let (name_b, comment_b) = match ws {
            Some(i) => (&header[..i], Some(&header[i + 1..])),
            None => (&header[..], None),
        };
        rec.name = String::from_utf8_lossy(name_b).into_owned();
        rec.comment = comment_b
            .filter(|c| !c.is_empty())
            .map(|c| String::from_utf8_lossy(c).into_owned());

        // 序列行：直到 '>'/'@'（FASTA/下一条）或 '+'（FASTQ）
        let mut sep = 0u8;
        let mut line = Vec::new();
        loop {
            match self.read_line_bytes(&mut line) {
                Ok(false) => {
                    self.done = true;
                    break;
                }
                Err(e) => return Some(Err(e)),
                Ok(true) => {}
            }
            if let Some(&first) = line.first() {
                if first == b'>' || first == b'@' {
                    self.pending_header = Some(line[1..].to_vec());
                    break;
                }
                if first == b'+' {
                    sep = b'+';
                    break;
                }
            }
            rec.seq.extend_from_slice(&line);
        }

        if sep == b'+' {
            // FASTQ：qual 读到与 seq 等长（kseq 允许跨行）
            while rec.qual.len() < rec.seq.len() {
                match self.read_line_bytes(&mut line) {
                    Ok(false) => {
                        self.done = true;
                        break;
                    }
                    Err(e) => return Some(Err(e)),
                    Ok(true) => {}
                }
                for &b in line.iter() {
                    if rec.qual.len() < rec.seq.len() {
                        rec.qual.push(b);
                    }
                }
            }
            if rec.qual.len() != rec.seq.len() {
                return Some(Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "qual 长度与 seq 不一致",
                )));
            }
        }
        Some(Ok(rec))
    }
}

impl<R: Read> Iterator for SeqIter<R> {
    type Item = std::io::Result<SeqRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        self.next_record()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_fasta_and_fastq() {
        let data = b">r1 a comment\nACGT\nACGT\n>r2\nTTTT\n@q1\nACGT\n+\nIIII\n@q2 long\nNNNN\n+\n!!!!\n";
        let mut it = SeqIter::new(&data[..]);
        let r1 = it.next().unwrap().unwrap();
        assert_eq!(r1.name, "r1");
        assert_eq!(r1.comment.as_deref(), Some("a comment"));
        assert_eq!(r1.seq, b"ACGTACGT");
        assert!(r1.qual.is_empty());
        let r2 = it.next().unwrap().unwrap();
        assert_eq!(r2.name, "r2");
        assert_eq!(r2.seq, b"TTTT");
        let q1 = it.next().unwrap().unwrap();
        assert_eq!(q1.name, "q1");
        assert_eq!(q1.seq, b"ACGT");
        assert_eq!(q1.qual, b"IIII");
        let q2 = it.next().unwrap().unwrap();
        assert_eq!(q2.name, "q2");
        assert_eq!(q2.comment.as_deref(), Some("long"));
        assert_eq!(q2.seq, b"NNNN");
        assert!(it.next().is_none());
    }
}
