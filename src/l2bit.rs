//! l2bit 参考基因组格式（`.l2b`），逐行对照 C 版 `l2bit.c`。
//!
//! 文件格式（小端）：
//! ```text
//! magic "L2B\1" | u32 dummy
//! u64 n_ctg, tot_len, n_ambi, n_mask, len_name, len_comm, n_pac
//! u64[n_ctg] 各 contig 长度（off 由累加得到）
//! (u64 st, u64 en)[n_ambi]  模糊碱基区间
//! (u64 st, u64 en)[n_mask]  软屏蔽（小写）区间
//! u64[n_pac]  2-bit 压缩序列，第 i 位 = pac[i>>5] >> ((i&31)*2) & 3
//! 名字串（NUL 分隔）| 注释串（NUL 分隔，空注释为单个 NUL）
//! ```

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;

use crate::dna::{splitmix64, NT4};

pub const MAGIC: &[u8; 4] = b"L2B\x01";

#[derive(Debug, Clone)]
pub struct Ctg {
    pub name: String,
    pub comm: Option<String>,
    pub len: u64,
    pub off: u64,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Intv {
    pub st: u64,
    pub en: u64,
}

/// 2-bit 参考基因组（对应 `l2b_t`）。
#[derive(Debug, Clone)]
pub struct L2b {
    pub tot_len: u64,
    pub ctg: Vec<Ctg>,
    pub ambi: Vec<Intv>,
    pub mask: Vec<Intv>,
    pub pac: Vec<u64>,
}

impl L2b {
    pub fn new() -> Self {
        L2b {
            tot_len: 0,
            ctg: Vec::new(),
            ambi: Vec::new(),
            mask: Vec::new(),
            pac: Vec::new(),
        }
    }

    /// 取全局坐标 i 处的 2-bit 碱基（对应 `l2b_get0`）。
    #[inline]
    pub fn get0(&self, i: u64) -> u8 {
        (self.pac[(i >> 5) as usize] >> ((i & 31) << 1)) as u8 & 3
    }

    /// C(1) → T(3)，对照 `l2b_c2t`。
    #[inline]
    pub(crate) fn c2t(b: u8) -> u8 {
        if b == 1 {
            3
        } else {
            b
        }
    }

    /// G(2) → A(0)，对照 `l2b_g2a`。
    #[inline]
    pub(crate) fn g2a(b: u8) -> u8 {
        if b == 2 {
            0
        } else {
            b
        }
    }

    /// 从 FASTA 导入（对应 `l2b_import`）。seed 用于模糊碱基的随机低 2 bit。
    pub fn import<R: Read>(reader: R, seed: u64) -> Result<Self> {
        let mut l2b = L2b::new();
        let mut rng = seed;
        for rec in FastaIter::new(reader) {
            let (name, comm, seq) = rec?;
            // l2b_format_seq：编码 + 模糊碱基随机化 + 小写标记
            let mut coded: Vec<u8> = Vec::with_capacity(seq.len());
            for &b in &seq {
                let mut c = NT4[b as usize];
                if c > 4 {
                    c = 4;
                }
                if c == 4 {
                    c |= (splitmix64(&mut rng) & 3) as u8;
                }
                if !b.is_ascii_uppercase() {
                    c |= 1 << 3;
                }
                coded.push(c);
            }
            l2b.add_seq(&coded, &name, comm.as_deref());
        }
        Ok(l2b)
    }

    /// 追加一条 contig（对应 `l2b_add_seq`）。
    fn add_seq(&mut self, seq: &[u8], name: &str, comm: Option<&str>) {
        let off = self.tot_len;
        let len = seq.len() as u64;
        self.ctg.push(Ctg {
            name: name.to_string(),
            comm: comm.map(str::to_string),
            len,
            off,
        });
        self.tot_len += len;

        let n_pac_old = self.pac.len();
        let n_pac_new = self.tot_len.div_ceil(32) as usize;
        self.pac.resize(n_pac_new, 0);
        let _ = n_pac_old; // C 版对新分配部分清零；Vec::resize 已保证

        let mut ambi_len = 0u64;
        let mut mask_len = 0u64;
        for (i, &c) in seq.iter().enumerate() {
            let x = off + i as u64;
            if c & 1 << 3 != 0 {
                mask_len += 1;
            } else if mask_len > 0 {
                self.mask.push(Intv {
                    st: x - mask_len,
                    en: x,
                });
                mask_len = 0;
            }
            if c & 1 << 2 != 0 {
                ambi_len += 1;
            } else if ambi_len > 0 {
                self.ambi.push(Intv {
                    st: x - ambi_len,
                    en: x,
                });
                ambi_len = 0;
            }
            self.pac[(x >> 5) as usize] |= ((c & 3) as u64) << ((x & 31) << 1);
        }
        // 注意：与 C 版一致，序列末尾未闭合的 ambi/mask 游程**不**落盘
    }

    /// 保存为 `.l2b`（对应 `l2b_save`）。
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path.as_ref()).with_context(|| format!("无法创建 {}", path.as_ref().display()))?;
        let mut w = BufWriter::new(file);
        let len_name: u64 = self.ctg.iter().map(|c| c.name.len() as u64 + 1).sum();
        let len_comm: u64 = self.ctg.iter().map(|c| c.comm.as_ref().map_or(1, |s| s.len() as u64 + 1)).sum();
        w.write_all(MAGIC)?;
        w.write_all(&0u32.to_le_bytes())?; // dummy
        w.write_all(&(self.ctg.len() as u64).to_le_bytes())?;
        w.write_all(&self.tot_len.to_le_bytes())?;
        w.write_all(&(self.ambi.len() as u64).to_le_bytes())?;
        w.write_all(&(self.mask.len() as u64).to_le_bytes())?;
        w.write_all(&len_name.to_le_bytes())?;
        w.write_all(&len_comm.to_le_bytes())?;
        w.write_all(&(self.pac.len() as u64).to_le_bytes())?;
        for c in &self.ctg {
            w.write_all(&c.len.to_le_bytes())?;
        }
        for iv in self.ambi.iter().chain(self.mask.iter()) {
            w.write_all(&iv.st.to_le_bytes())?;
            w.write_all(&iv.en.to_le_bytes())?;
        }
        for &x in &self.pac {
            w.write_all(&x.to_le_bytes())?;
        }
        for c in &self.ctg {
            w.write_all(c.name.as_bytes())?;
            w.write_all(&[0])?;
        }
        for c in &self.ctg {
            if let Some(comm) = &c.comm {
                w.write_all(comm.as_bytes())?;
            }
            w.write_all(&[0])?;
        }
        w.flush()?;
        Ok(())
    }

    /// 从 `.l2b` 加载（对应 `l2b_load`）。
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let data = std::fs::read(path.as_ref()).with_context(|| format!("无法读取 {}", path.as_ref().display()))?;
        Self::from_bytes(&data).with_context(|| format!("解析 .l2b 失败: {}", path.as_ref().display()))
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut r = Cursor::new(data);
        let magic = r.take(4)?;
        if magic != MAGIC {
            bail!("magic 不匹配（期望 L2B\\1）");
        }
        let _dummy = r.u32()?;
        let n_ctg = r.u64()?;
        let tot_len = r.u64()?;
        let n_ambi = r.u64()?;
        let n_mask = r.u64()?;
        let len_name = r.u64()?;
        let len_comm = r.u64()?;
        let n_pac = r.u64()?;

        // 头部计数全部是文件可控值。分配前按固定记录尺寸把它们与文件
        // 剩余字节核对，损坏的 .l2b 不能靠声明值驱动巨量分配（旧路径
        // Vec::with_capacity(1<<60) 直接 capacity overflow panic）。
        let remaining = r.rest().len() as u64;
        let needed = n_ctg
            .checked_mul(8)
            .and_then(|v| v.checked_add(n_ambi.checked_mul(16)?))
            .and_then(|v| v.checked_add(n_mask.checked_mul(16)?))
            .and_then(|v| v.checked_add(n_pac.checked_mul(8)?))
            .and_then(|v| v.checked_add(len_name))
            .and_then(|v| v.checked_add(len_comm));
        if needed.is_none_or(|n| n > remaining) {
            bail!("头部尺寸与文件剩余字节不一致（需要 {needed:?}，剩余 {remaining}）");
        }
        let n_ctg = usize::try_from(n_ctg).context("n_ctg 超出平台大小")?;
        let n_ambi = usize::try_from(n_ambi).context("n_ambi 超出平台大小")?;
        let n_mask = usize::try_from(n_mask).context("n_mask 超出平台大小")?;
        let len_name = usize::try_from(len_name).context("len_name 超出平台大小")?;
        let len_comm = usize::try_from(len_comm).context("len_comm 超出平台大小")?;
        let n_pac = usize::try_from(n_pac).context("n_pac 超出平台大小")?;

        // 压缩序列长度由 tot_len 唯一决定；不一致会让 get0 越界或读错位置。
        if n_pac as u64 != tot_len.div_ceil(32) {
            bail!("n_pac={n_pac} 与 tot_len={tot_len} 不一致");
        }

        let mut ctg = Vec::with_capacity(n_ctg);
        let mut off = 0u64;
        for _ in 0..n_ctg {
            let len = r.u64()?;
            ctg.push(Ctg {
                name: String::new(),
                comm: None,
                len,
                off,
            });
            off = off.checked_add(len).context("contig 长度累加溢出")?;
        }
        if off != tot_len {
            bail!("contig 长度之和与 tot_len 不一致");
        }

        let read_intvs = |r: &mut Cursor, n: usize| -> Result<Vec<Intv>> {
            let mut v = Vec::with_capacity(n);
            for _ in 0..n {
                let st = r.u64()?;
                let en = r.u64()?;
                v.push(Intv { st, en });
            }
            Ok(v)
        };
        let ambi = read_intvs(&mut r, n_ambi)?;
        let mask = read_intvs(&mut r, n_mask)?;

        let mut pac = vec![0u64; n_pac];
        for x in pac.iter_mut() {
            *x = r.u64()?;
        }
        let names = r.take(len_name)?;
        let comms = r.take(len_comm)?;
        if !r.rest().is_empty() {
            bail!("文件末尾有多余数据");
        }

        // 同步名字与注释（对应 l2b_load 末尾循环）
        let mut p_name = 0usize;
        let mut p_comm = 0usize;
        for c in ctg.iter_mut() {
            let end = names[p_name..].iter().position(|&b| b == 0).context("名字串未以 NUL 结尾")?;
            c.name = String::from_utf8(names[p_name..p_name + end].to_vec())?;
            p_name += end + 1;
            if comms[p_comm] != 0 {
                let cend = comms[p_comm..].iter().position(|&b| b == 0).context("注释串未以 NUL 结尾")?;
                c.comm = Some(String::from_utf8(comms[p_comm..p_comm + cend].to_vec())?);
                p_comm += cend + 1;
            } else {
                p_comm += 1;
            }
        }
        if p_name != len_name || p_comm != len_comm {
            bail!("名字/注释串长度与头部不一致");
        }

        Ok(L2b {
            tot_len,
            ctg,
            ambi,
            mask,
            pac,
        })
    }

    /// 区间 [st, en) 反查 contig（对应 `l2b_intv2cid`）。
    /// 返回 contig id 与 contig 内起点；rev 表示是否落在反向互补拷贝上。
    /// 错误码与 C 版一致：-2 跨 contig 边界，-3 超出总长。
    pub fn intv2cid(&self, st: u64, en: u64, rev: &mut bool) -> (i64, i64) {
        assert!(st < en);
        if en > self.tot_len * 2 {
            return (-3, 0);
        }
        if st < self.tot_len && self.tot_len < en {
            return (-2, 0);
        }
        *rev = st >= self.tot_len;
        let s = if st < self.tot_len { st } else { self.tot_len * 2 - en };
        self.pos2cid(s, en - st)
    }

    fn pos2cid(&self, s: u64, len: u64) -> (i64, i64) {
        let mut lo = 0usize;
        let mut hi = self.ctg.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            let c = &self.ctg[mid];
            if c.off <= s && s < c.off + c.len {
                let cst = (s - c.off) as i64;
                let id = if s + len <= c.off + c.len { mid as i64 } else { -1 };
                return (id, cst);
            } else if s < c.off {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        (-1, 0)
    }

    /// 取 contig tid 的 [st, en) 序列（2-bit，模糊碱基置 4），对应 `l2b_getseq`。
    pub fn getseq(&self, tid: i64, mut st: i64, mut en: i64, seq: &mut [u8]) -> i64 {
        if tid < 0 || tid as usize >= self.ctg.len() {
            return -1;
        }
        let ctg = &self.ctg[tid as usize];
        if st < 0 {
            st = 0;
        }
        if en > ctg.len as i64 {
            en = ctg.len as i64;
        }
        let gst = ctg.off as i64 + st;
        let gen = ctg.off as i64 + en;
        for (i, g) in (gst..gen).enumerate() {
            seq[i] = self.get0(g as u64);
        }
        // 模糊碱基区间覆盖为 4
        let (aid, n_ambi) = self.getambi(tid, st, en);
        if aid >= 0 {
            for iv in &self.ambi[aid as usize..aid as usize + n_ambi as usize] {
                let mut s = iv.st as i64;
                let mut e = iv.en as i64;
                if s < gst {
                    s = gst;
                }
                if e > gen {
                    e = gen;
                }
                for x in s..e {
                    seq[(x - gst) as usize] = 4;
                }
            }
        }
        gen - gst
    }

    /// 查询与 [st, en) 相交的 ambi 区间（对应 `l2b_getambi`），返回 (起始下标, 个数)。
    pub fn getambi(&self, tid: i64, mut st: i64, mut en: i64) -> (i64, i32) {
        if tid < 0 || tid as usize >= self.ctg.len() {
            return (-1, 0);
        }
        let ctg = &self.ctg[tid as usize];
        if st < 0 {
            st = 0;
        }
        if en > ctg.len as i64 {
            en = ctg.len as i64;
        }
        if st >= en {
            return (-1, 0);
        }
        let g_beg = ctg.off as i64 + st;
        let g_end = ctg.off as i64 + en;

        let mut lo = 0usize;
        let mut hi = self.ambi.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.ambi[mid].en as i64 > g_beg {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        let i_st = lo;

        hi = self.ambi.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.ambi[mid].st as i64 >= g_end {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        let i_en = lo;
        if i_en == i_st {
            return (-1, 0);
        }
        (i_st as i64, (i_en - i_st) as i32)
    }
}

impl Default for L2b {
    fn default() -> Self {
        Self::new()
    }
}

/// 从 `.l2b` 还原参考序列并打印（对应 C 版 `main_getref`）。
pub fn getref<P: AsRef<Path>>(path: P) -> Result<()> {
    let l2b = L2b::load(path)?;
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    const BASES: &[u8; 5] = b"ACGTN";
    for (i, ctg) in l2b.ctg.iter().enumerate() {
        writeln!(out, ">{}", ctg.name)?;
        let mut seq = vec![0u8; ctg.len as usize];
        l2b.getseq(i as i64, 0, ctg.len as i64, &mut seq);
        for b in seq.iter_mut() {
            *b = BASES[*b as usize];
        }
        out.write_all(&seq)?;
        writeln!(out)?;
    }
    out.flush()?;
    Ok(())
}

pub struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        Cursor { data, pos: 0 }
    }

    pub(crate) fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            bail!("文件截断");
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    pub(crate) fn u32(&mut self) -> Result<u32> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes(b.try_into().unwrap()))
    }

    pub(crate) fn u64(&mut self) -> Result<u64> {
        let b = self.take(8)?;
        Ok(u64::from_le_bytes(b.try_into().unwrap()))
    }

    pub(crate) fn rest(&self) -> &'a [u8] {
        &self.data[self.pos..]
    }
}

/// kseq 兼容的极简 FASTA 迭代器：
/// name 取到第一个空白字符；comment 为紧随其后的行内剩余部分（前导空白保留）；
/// 序列行只跳过 '\n'（与 C 版 kseq 一致，'\r' 会被保留）。
struct FastaIter<R: Read> {
    reader: BufReader<R>,
    pending_header: Option<Vec<u8>>,
    done: bool,
}

impl<R: Read> FastaIter<R> {
    fn new(reader: R) -> Self {
        FastaIter {
            reader: BufReader::new(reader),
            pending_header: None,
            done: false,
        }
    }
}

impl<R: Read> Iterator for FastaIter<R> {
    /// (name, comment, seq)。seq 为原始 ASCII 字节。
    type Item = Result<(String, Option<String>, Vec<u8>)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        // 找到下一条头行（'>' 或 '@'，与 kseq 一致）
        let header = if let Some(h) = self.pending_header.take() {
            h
        } else {
            let mut line = Vec::new();
            loop {
                line.clear();
                match self.reader.read_until(b'\n', &mut line) {
                    Ok(0) => {
                        self.done = true;
                        return None;
                    }
                    Err(e) => return Some(Err(e.into())),
                    Ok(_) => {}
                }
                if line.first() == Some(&b'>') || line.first() == Some(&b'@') {
                    break;
                }
            }
            line
        };

        // 解析 name / comment（kseq：name 到第一个空白；comment 为行内剩余）
        let body = &header[1..];
        let body = match body.iter().position(|&b| b == b'\n') {
            Some(i) => &body[..i],
            None => body,
        };
        let ws = body.iter().position(|&b| matches!(b, b' ' | b'\t' | b'\r' | 0x0b | 0x0c));
        let (name_bytes, comm_bytes) = match ws {
            Some(i) => (&body[..i], Some(&body[i + 1..])),
            None => (body, None),
        };
        let name = String::from_utf8_lossy(name_bytes).into_owned();
        let comm = comm_bytes.map(|b| String::from_utf8_lossy(b).into_owned());

        // 读序列行直到下一个头
        let mut seq = Vec::new();
        let mut line = Vec::new();
        loop {
            line.clear();
            match self.reader.read_until(b'\n', &mut line) {
                Ok(0) => {
                    self.done = true;
                    break;
                }
                Err(e) => return Some(Err(e.into())),
                Ok(_) => {}
            }
            if line.first() == Some(&b'>') || line.first() == Some(&b'@') {
                self.pending_header = Some(line.clone());
                break;
            }
            // kseq 语义：跳过 '\n'，其余字节（含 '\r'）全部保留
            for &b in line.iter() {
                if b == b'\n' {
                    break;
                }
                seq.push(b);
            }
        }
        Some(Ok((name, comm, seq)))
    }
}

/// 打开 FASTA（透明支持 gzip，对应 C 版 gzopen）。
pub fn open_fasta(path: &str) -> Result<Box<dyn Read>> {
    if path == "-" {
        return Ok(Box::new(std::io::stdin()));
    }
    let file = File::open(path).with_context(|| format!("无法打开 {path}"))?;
    let mut br = BufReader::new(file);
    // fill_buf 只peek不消费，避免吃掉序列头字节
    let is_gz = {
        let buf = br.fill_buf().context("读取文件头失败")?;
        buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
    };
    if is_gz {
        Ok(Box::new(MultiGzDecoder::new(br)))
    } else {
        Ok(Box::new(br))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn import_encodes_bases_and_runs() {
        // 末尾补一个 A，保证小写游程不在序列末尾（末尾游程不落盘，见下一测试）
        let fa = ">c1\nACGTNnacgtA\n";
        let l2b = L2b::import(fa.as_bytes(), 11).unwrap();
        assert_eq!(l2b.tot_len, 11);
        assert_eq!(l2b.ctg.len(), 1);
        // N 与 n 都是模糊碱基：区间 [4,6)（大写 N 在 4，小写 n 在 5）
        assert_eq!(l2b.ambi.len(), 1);
        assert_eq!((l2b.ambi[0].st, l2b.ambi[0].en), (4, 6));
        // 小写 nacgt 区间 [5,10) 为软屏蔽
        assert_eq!(l2b.mask.len(), 1);
        assert_eq!((l2b.mask[0].st, l2b.mask[0].en), (5, 10));
        // 确定碱基编码正确
        assert_eq!(l2b.get0(0), 0); // A
        assert_eq!(l2b.get0(1), 1); // C
        assert_eq!(l2b.get0(2), 2); // G
        assert_eq!(l2b.get0(3), 3); // T
        assert_eq!(l2b.get0(6), 0); // a
    }

    #[test]
    fn trailing_ambi_run_not_recorded() {
        // 与 C 版行为一致：序列末尾未闭合的模糊游程不落盘
        let fa = ">c1\nACGTNNN\n";
        let l2b = L2b::import(fa.as_bytes(), 11).unwrap();
        assert!(l2b.ambi.is_empty());
    }

    #[test]
    fn save_load_roundtrip() {
        let fa = ">c1 some comment\nACGTNNacgt\n>c2\nTTTT\n";
        let l2b = L2b::import(fa.as_bytes(), 11).unwrap();
        let bytes = {
            let tmp = std::env::temp_dir().join("minibwa-rust-test.l2b");
            l2b.save(&tmp).unwrap();
            std::fs::read(&tmp).unwrap()
        };
        let l2b2 = L2b::from_bytes(&bytes).unwrap();
        assert_eq!(l2b2.tot_len, l2b.tot_len);
        assert_eq!(l2b2.ctg.len(), 2);
        assert_eq!(l2b2.ctg[0].name, "c1");
        assert_eq!(l2b2.ctg[0].comm.as_deref(), Some("some comment"));
        assert_eq!(l2b2.ctg[1].name, "c2");
        assert_eq!(l2b2.ctg[1].comm, None);
        assert_eq!(l2b2.pac, l2b.pac);
        assert_eq!(l2b2.ambi.len(), l2b.ambi.len());
    }

    #[test]
    fn getseq_masks_ambi() {
        let fa = ">c1\nACGTNNACGT\n";
        let l2b = L2b::import(fa.as_bytes(), 11).unwrap();
        let mut seq = vec![0u8; 10];
        l2b.getseq(0, 0, 10, &mut seq);
        assert_eq!(&seq[..4], &[0, 1, 2, 3]);
        assert_eq!(&seq[4..6], &[4, 4]);
        assert_eq!(&seq[6..], &[0, 1, 2, 3]);
    }

    #[test]
    fn intv2cid_forward_and_rev() {
        let fa = ">c1\nACGTACGT\n>c2\nTTTT\n";
        let l2b = L2b::import(fa.as_bytes(), 11).unwrap();
        let mut rev = false;
        let (id, cst) = l2b.intv2cid(2, 5, &mut rev);
        assert_eq!((id, cst, rev), (0, 2, false));
        // 反向拷贝：[tot+4, tot+8) 映射回正向 [4,8)，落在 c1
        let (id, cst) = l2b.intv2cid(l2b.tot_len + 4, l2b.tot_len + 8, &mut rev);
        assert_eq!((id, cst, rev), (0, 4, true));
        // 跨 contig 边界：pos2cid 返回 -1
        let (id, _) = l2b.intv2cid(6, 10, &mut rev);
        assert_eq!(id, -1);
        // 跨越正反向连接处（tot_len）：返回 -2
        let (id, _) = l2b.intv2cid(10, 14, &mut rev);
        assert_eq!(id, -2);
    }
}
