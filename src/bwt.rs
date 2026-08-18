//! BWT + 采样后缀数组，逐行对照 C 版 `bwt.c`。
//!
//! `.mbw` 文件格式（小端）：
//! ```text
//! magic "MBW\2" | u32 sa_bit
//! u64 primary | u64[4] L2[1..5]
//! u64[data_len] BWT 数据块（每块 64 字节：u64[4] 计数 + u32[8] 128nt 2-bit）
//! u64 n_sa | u64[n_sa] 采样 SA
//! ```

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

pub const MAGIC: &[u8; 4] = b"MBW\x02";

const CNT_SHIFT: u32 = 56;
const CNT_MASK: u64 = (1u64 << CNT_SHIFT) - 1;

/// SA 区间（对应 `mb_sai_t`）。x[0]/x[1] 为前后向搜索坐标，半开半合。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Sai {
    pub x: [u64; 2],
    pub size: u64,
    pub info: u64,
}

/// BWT 索引（对应 `mb_bwt_t`）。
pub struct Bwt {
    pub primary: u64,
    pub l2: [u64; 5],
    pub seq_len: u64,
    pub data: Vec<u64>,
    cnt_table: [u32; 256],
    pub pre_len: u32,
    pub pre: Option<Vec<Sai>>,
    pub sa_bit: u32,
    pub sa: Vec<u64>,
}

fn gen_cnt_table() -> [u32; 256] {
    let mut cnt = [0u32; 256];
    for (i, x) in cnt.iter_mut().enumerate() {
        let mut v = 0u32;
        for j in 0..4u32 {
            let c = ((i & 3) == j as usize) as u32
                + (((i >> 2) & 3) == j as usize) as u32
                + (((i >> 4) & 3) == j as usize) as u32
                + ((i >> 6) == j as usize) as u32;
            v |= c << (j << 3);
        }
        *x = v;
    }
    cnt
}

/// data 长度（u64 个数），对应 `mb_bwt_data_len`。
pub fn data_len(len: u64) -> u64 {
    let bwt_len = len.div_ceil(128) * 4;
    let occ_len = (len.div_ceil(128) + 1) * 4; // +1 存末尾总计数
    bwt_len + occ_len
}

impl Bwt {
    pub fn new() -> Self {
        Bwt {
            primary: 0,
            l2: [0; 5],
            seq_len: 0,
            data: Vec::new(),
            cnt_table: gen_cnt_table(),
            pre_len: 0,
            pre: None,
            sa_bit: u32::MAX,
            sa: Vec::new(),
        }
    }

    /// 从 2-bit 编码的原始 BWT 串构建块状布局（对应 `mb_bwt_init_from_raw`）。
    /// raw[i] ∈ 0..3，长度 len；primary = S^{-1}(0)。
    pub fn from_raw(raw: &[u8], len: u64, primary: u64) -> Self {
        let mut bwt = Bwt::new();
        bwt.primary = primary;
        bwt.seq_len = len;
        let dl = data_len(len) as usize;
        bwt.data = vec![0u64; dl];
        if len == 0 {
            return bwt;
        }

        let mut c = [0u64; 4];
        let mut x = [0u64; 4];
        let mut k = 0usize;
        let mut last_c = usize::MAX;
        for i in 0..len as usize {
            let a = raw[i] & 3;
            if i & 0x7f == 0 {
                // 每块 128bp：先落上一块的 128nt 数据，再记本块起始计数
                if i > 0 {
                    bwt.data[k..k + 4].copy_from_slice(&x);
                    k += 4;
                }
                last_c = k;
                bwt.data[k..k + 4].copy_from_slice(&c);
                k += 4;
                x = [0u64; 4];
            } else if i & 0x3f == 0 {
                // 块内第 64nt：把前 64nt 的计数差写入高 8 bit
                for (slot, &cj) in bwt.data[last_c..last_c + 4].iter_mut().zip(&c) {
                    *slot |= (cj - *slot) << CNT_SHIFT;
                }
            }
            c[a as usize] += 1;
            x[(i & 0x7f) >> 5] |= (a as u64) << ((i & 0x1f) << 1);
        }
        // 末块：数据 + 总计数
        bwt.data[k..k + 4].copy_from_slice(&x);
        k += 4;
        bwt.data[k..k + 4].copy_from_slice(&c);
        k += 4;
        assert_eq!(k, dl);
        bwt.l2[0] = 0;
        for (i, &ci) in c.iter().enumerate() {
            bwt.l2[i + 1] = bwt.l2[i] + ci;
        }
        assert_eq!(bwt.l2[4], len);
        bwt
    }

    /// 取 $-removed BWT 串第 k 位碱基（对应 `bwt_B0`）。
    #[inline]
    fn b0(&self, k: u64) -> u8 {
        ((self.data[((k >> 7) << 3) as usize + 4 + (((k & 127) >> 5) as usize)] >> ((k & 31) << 1)) & 3) as u8
    }

    #[inline]
    fn rank_aux1(y: u64, c: u8) -> u64 {
        // 将核苷酸计数归约为位计数（对应 rank_aux1）
        let y = (((if c & 2 != 0 { y } else { !y }) >> 1) & (if c & 1 != 0 { y } else { !y })) & 0x5555_5555_5555_5555;
        y.count_ones() as u64
    }

    #[inline]
    fn rank_aux4(&self, x: u32) -> u32 {
        self.cnt_table[(x & 0xff) as usize]
            + self.cnt_table[((x >> 8) & 0xff) as usize]
            + self.cnt_table[((x >> 16) & 0xff) as usize]
            + self.cnt_table[(x >> 24) as usize]
    }

    /// 单碱基 rank：B[0..k) 中碱基 c 的个数（对应 `mb_bwt_rank11`）。
    pub fn rank11(&self, mut k: u64, c: u8) -> u64 {
        if k == 0 {
            return 0;
        }
        if k == self.seq_len + 1 {
            return self.l2[c as usize + 1] - self.l2[c as usize];
        }
        k -= 1;
        if k >= self.primary {
            k -= 1; // $ 不在 BWT 串中
        }
        let mask = if k & 0x7f >= 64 {
            (1u64 << (64 - CNT_SHIFT)) - 1
        } else {
            0
        };
        let block = ((k >> 7) << 3) as usize;
        let mut n = (self.data[block + c as usize] & CNT_MASK) + ((self.data[block + c as usize] >> CNT_SHIFT) & mask);
        let mut p = block + 4 + if k & 0x7f >= 64 { 2 } else { 0 };
        let end = block + 4 + ((k & 0x7f) >> 5) as usize;
        if p < end {
            n += Self::rank_aux1(self.data[p], c);
            p += 1;
        }
        n += Self::rank_aux1(self.data[p] << ((!k & 0x1f) << 1), c);
        if c == 0 {
            n -= !k & 0x1f; // 左移补零会被误计为 A，此处修正
        }
        n
    }

    /// 四碱基 rank：B[0..k) 中各碱基个数（对应 `mb_bwt_rank1a`）。
    pub fn rank1a(&self, mut k: u64, cnt: &mut [u64; 4]) {
        if k == 0 {
            *cnt = [0; 4];
            return;
        }
        k -= 1;
        if k >= self.primary {
            k -= 1;
        }
        let block = ((k >> 7) << 3) as usize;
        let mask = if k & 0x7f >= 64 {
            (1u64 << (64 - CNT_SHIFT)) - 1
        } else {
            0
        };
        for (j, cnt_j) in cnt.iter_mut().enumerate() {
            *cnt_j = (self.data[block + j] & CNT_MASK) + ((self.data[block + j] >> CNT_SHIFT) & mask);
        }
        // 块内 2-bit 数据按 u32 计（小端）
        let data32 = &self.data[block + 4..];
        let data32: Vec<u32> = data32.iter().flat_map(|&w| [w as u32, (w >> 32) as u32]).collect();
        let mut q = if k & 0x7f >= 64 { 4 } else { 0 };
        let end = ((k & 0x7f) >> 4) as usize;
        let mut x = 0u32;
        while q < end {
            x = x.wrapping_add(self.rank_aux4(data32[q]));
            q += 1;
        }
        let tmp = data32[q] << ((!k & 0xf) << 1);
        x = x.wrapping_add(self.rank_aux4(tmp).wrapping_sub((!k & 0xf) as u32));
        cnt[0] += (x & 0xff) as u64;
        cnt[1] += ((x >> 8) & 0xff) as u64;
        cnt[2] += ((x >> 16) & 0xff) as u64;
        cnt[3] += (x >> 24) as u64;
    }

    /// 双坐标四碱基 rank（对应 `mb_bwt_rank2a`）。
    pub fn rank2a(&self, k: u64, l: u64, cntk: &mut [u64; 4], cntl: &mut [u64; 4]) {
        let mut k1 = k - 1;
        let mut l1 = l - 1;
        if k1 >= self.primary {
            k1 -= 1;
        }
        if l1 >= self.primary {
            l1 -= 1;
        }
        if k1 >> 7 != l1 >> 7 || k == 0 || l == 0 {
            self.rank1a(k, cntk);
            self.rank1a(l, cntl);
        } else if l - k == 1 {
            let z = k - u64::from(k > self.primary);
            self.rank1a(k, cntk);
            *cntl = *cntk;
            cntl[self.b0(z) as usize] += 1;
        } else {
            let block = ((k1 >> 7) << 3) as usize;
            let mask = if k1 & 0x7f >= 64 {
                (1u64 << (64 - CNT_SHIFT)) - 1
            } else {
                0
            };
            for (j, cntk_j) in cntk.iter_mut().enumerate() {
                *cntk_j = (self.data[block + j] & CNT_MASK) + ((self.data[block + j] >> CNT_SHIFT) & mask);
            }
            let data32: Vec<u32> = self.data[block + 4..block + 8]
                .iter()
                .flat_map(|&w| [w as u32, (w >> 32) as u32])
                .collect();
            let mut q = if k1 & 0x7f >= 64 { 4 } else { 0 };
            let endk = ((k1 & 0x7f) >> 4) as usize;
            let endl = ((l1 & 0x7f) >> 4) as usize;
            let mut x = 0u32;
            while q < endk {
                x = x.wrapping_add(self.rank_aux4(data32[q]));
                q += 1;
            }
            let mut y = x;
            let tmp = data32[q] << ((!k1 & 0xf) << 1);
            x = x.wrapping_add(self.rank_aux4(tmp).wrapping_sub((!k1 & 0xf) as u32));
            while q < endl {
                y = y.wrapping_add(self.rank_aux4(data32[q]));
                q += 1;
            }
            let tmp = data32[q] << ((!l1 & 0xf) << 1);
            y = y.wrapping_add(self.rank_aux4(tmp).wrapping_sub((!l1 & 0xf) as u32));
            *cntl = *cntk;
            cntk[0] += (x & 0xff) as u64;
            cntk[1] += ((x >> 8) & 0xff) as u64;
            cntk[2] += ((x >> 16) & 0xff) as u64;
            cntk[3] += (x >> 24) as u64;
            cntl[0] += (y & 0xff) as u64;
            cntl[1] += ((y >> 8) & 0xff) as u64;
            cntl[2] += ((y >> 16) & 0xff) as u64;
            cntl[3] += (y >> 24) as u64;
        }
    }

    /// 双向 BWT 扩展（对应 `mb_bwt_extend`）。is_back=true 为向后（5'→3' 反向搜索）。
    pub fn extend(&self, ik: &Sai, ok: &mut [Sai; 4], is_back: bool) {
        let (d, nd) = if is_back { (1, 0) } else { (0, 1) };
        let mut tk = [0u64; 4];
        let mut tl = [0u64; 4];
        self.rank2a(ik.x[nd], ik.x[nd] + ik.size, &mut tk, &mut tl);
        for i in 0..4 {
            ok[i].x[nd] = self.l2[i] + 1 + tk[i]; // +1 补偿被移除的 $
            ok[i].size = tl[i] - tk[i];
            tl[i] -= tk[i];
        }
        let cross = u64::from(ik.x[nd] <= self.primary && ik.x[nd] + ik.size > self.primary);
        ok[3].x[d] = ik.x[d] + cross;
        ok[2].x[d] = ok[3].x[d] + tl[3];
        ok[1].x[d] = ok[2].x[d] + tl[2];
        ok[0].x[d] = ok[1].x[d] + tl[1];
    }

    /// 单碱基初始区间（对应 `mb_bwt_set_intv`）。
    pub fn set_intv(&self, c: u8) -> Sai {
        Sai {
            x: [self.l2[c as usize] + 1, self.l2[(3 - c) as usize] + 1],
            size: self.l2[c as usize + 1] - self.l2[c as usize],
            info: 0,
        }
    }

    /// 逆 CSA 一步（对应 `bwt_invPsi`）。
    #[inline]
    fn inv_psi(&self, k: u64) -> u64 {
        let x = k - u64::from(k > self.primary);
        let c = self.b0(x);
        if k == self.primary {
            return 0;
        }
        self.l2[c as usize] + 1 + self.rank11(k, c)
    }

    /// 生成采样 SA（对应 `mb_bwt_gen_sa`）。
    pub fn gen_sa(&mut self, sa_bit: u32) {
        assert!(!self.data.is_empty());
        self.sa_bit = sa_bit;
        let n_sa = ((self.seq_len + (1 << sa_bit)) >> sa_bit) as usize;
        self.sa = vec![0u64; n_sa];
        let mask = (1u64 << sa_bit) - 1;
        let mut isa = 0u64;
        let mut sa = self.seq_len;
        for _ in 0..self.seq_len {
            if isa & mask == 0 {
                self.sa[(isa >> sa_bit) as usize] = sa;
            }
            sa -= 1;
            isa = self.inv_psi(isa);
        }
        if isa & mask == 0 {
            self.sa[(isa >> sa_bit) as usize] = sa;
        }
        self.sa[0] = u64::MAX;
    }

    /// 由 BWT 行号求 SA 值（对应 `mb_bwt_sa`）。
    pub fn sa(&self, mut k: u64) -> u64 {
        let mut steps = 0u64;
        let mask = (1u64 << self.sa_bit) - 1;
        while k & mask != 0 {
            steps += 1;
            k = self.inv_psi(k);
        }
        // C 版依赖无符号回绕：sa[0] 被置为 -1，走经 0 行时 steps + (-1) 恰好
        // 等于正确 SA 值 mod (seq_len+1)。Rust 用 wrapping_add 复现该语义。
        steps.wrapping_add(self.sa[(k >> self.sa_bit) as usize])
    }

    /// 保存为 `.mbw`（对应 `mb_bwt_save`）。
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path.as_ref()).with_context(|| format!("无法创建 {}", path.as_ref().display()))?;
        let mut w = BufWriter::new(file);
        w.write_all(MAGIC)?;
        w.write_all(&self.sa_bit.to_le_bytes())?;
        w.write_all(&self.primary.to_le_bytes())?;
        for &v in &self.l2[1..] {
            w.write_all(&v.to_le_bytes())?;
        }
        for &v in &self.data {
            w.write_all(&v.to_le_bytes())?;
        }
        let n_sa = self.sa.len() as u64;
        w.write_all(&n_sa.to_le_bytes())?;
        if self.sa_bit != u32::MAX && n_sa > 0 {
            for &v in &self.sa {
                w.write_all(&v.to_le_bytes())?;
            }
        }
        w.flush()?;
        Ok(())
    }

    /// 从 `.mbw` 加载（对应 `mb_bwt_load`）。
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let data = std::fs::read(path.as_ref()).with_context(|| format!("无法读取 {}", path.as_ref().display()))?;
        let mut r = super::l2bit::Cursor::new(&data);
        if r.take(4)? != MAGIC {
            bail!("magic 不匹配（期望 MBW\\2）: {}", path.as_ref().display());
        }
        let mut bwt = Bwt::new();
        bwt.sa_bit = {
            let b = r.take(4)?;
            u32::from_le_bytes(b.try_into().unwrap())
        };
        let mut x = [0u64; 5];
        for v in x.iter_mut() {
            let b = r.take(8)?;
            *v = u64::from_le_bytes(b.try_into().unwrap());
        }
        bwt.primary = x[0];
        bwt.l2[0] = 0;
        bwt.l2[1..].copy_from_slice(&x[1..]);
        bwt.seq_len = bwt.l2[4];
        let dl = data_len(bwt.seq_len) as usize;
        bwt.data = vec![0u64; dl];
        for v in bwt.data.iter_mut() {
            let b = r.take(8)?;
            *v = u64::from_le_bytes(b.try_into().unwrap());
        }
        let b = r.take(8)?;
        let n_sa = u64::from_le_bytes(b.try_into().unwrap());
        if bwt.sa_bit != u32::MAX && n_sa > 0 {
            let expected = (bwt.seq_len + (1u64 << bwt.sa_bit)) >> bwt.sa_bit;
            if n_sa != expected {
                bail!("n_sa 与 seq_len/sa_bit 不一致");
            }
            bwt.sa = vec![0u64; n_sa as usize];
            for v in bwt.sa.iter_mut() {
                let b = r.take(8)?;
                *v = u64::from_le_bytes(b.try_into().unwrap());
            }
        }
        if !r.rest().is_empty() {
            bail!("文件末尾有多余数据");
        }
        Ok(bwt)
    }

    /************** SMEM（对应 mb_bwt_back / mb_bwt_smem） **************/

    /// 从 pos 向后搜索（对应 `mb_bwt_back`）。返回停止位置 i。
    fn back(&self, q: &[u8], st: i64, pos: i64, min_occ: i64, p: &mut Sai) -> i64 {
        let mut i = pos - 1;
        if self.pre.is_some() && pos - st >= self.pre_len as i64 {
            // 用预算 kmer 缓存代替逐碱基扩展
            let mut z = 0usize;
            let mut l = 0i64;
            i = pos;
            while l < self.pre_len as i64 {
                z = z << 2 | q[i as usize] as usize;
                i -= 1;
                l += 1;
            }
            // 循环结束时 i == pos - pre_len，与 C 版 `for (i = pos; l < pre_len; --i)` 一致
            *p = self.pre.as_ref().unwrap()[z];
        } else {
            p.size = 0;
        }
        if p.size < min_occ as u64 {
            // kmer 跳得太远或无缓存：回到标准流程
            *p = self.set_intv(q[pos as usize]);
            i = pos - 1;
        }
        let mut ok = [Sai::default(); 4];
        while i >= st {
            let c = q[i as usize];
            if c > 3 {
                break;
            }
            self.extend(p, &mut ok, true);
            if ok[c as usize].size < min_occ as u64 {
                break;
            }
            *p = ok[c as usize];
            i -= 1;
        }
        i
    }

    /// 查找 super MEM（对应 `mb_bwt_smem`）。
    /// q 为 2-bit 编码查询序列；从 x 开始找长度 ≥ min_len 的 SMEM。
    /// 结果写入 p（info = x<<32 | end），返回下一轮起点。
    pub fn smem(&self, q: &[u8], x: i64, min_len: i64, min_occ: i64, p: &mut Sai) -> i64 {
        let len = q.len() as i64;
        p.size = 0;
        if len - x < min_len {
            return len;
        }
        let mut xn = -1i64;
        for i in x..x + min_len {
            if q[i as usize] > 3 {
                xn = i;
            }
        }
        if xn >= 0 {
            return xn + 1;
        }
        let mut ik = Sai::default();
        let mut i = self.back(q, x, x + min_len - 1, min_occ, &mut ik);
        if i >= x {
            return i + 1;
        }
        let mut j = x + min_len;
        let mut ok = [Sai::default(); 4];
        while j < len {
            if q[j as usize] > 3 {
                break;
            }
            let c = 3 - q[j as usize];
            self.extend(&ik, &mut ok, false);
            if ok[c as usize].size < min_occ as u64 {
                break;
            }
            ik = ok[c as usize];
            j += 1;
        }
        *p = ik;
        p.info = (x as u64) << 32 | j as u64;
        if j == len {
            return len;
        }
        i = if q[j as usize] > 3 {
            j
        } else {
            self.back(q, x + 1, j, min_occ, &mut ik)
        };
        i + 1
    }

    /// 构建 kmer 区间缓存（对应 `mb_bwt_cache`），len ≤ 15。
    pub fn cache(&mut self, len: u32) {
        assert!(len <= 15);
        self.pre_len = len;
        let mut pre = vec![Sai::default(); 1 << (len * 2)];
        self.count_kmer(len, &mut pre);
        self.pre = Some(pre);
    }

    /// kmer 计数（对应 `mb_bwt_count_kmer`，改编自 ropebwt3 kount）。
    fn count_kmer(&self, depth: u32, s: &mut [Sai]) {
        struct Frame {
            p: Sai,
            d: u32,
            c: u8,
        }
        let mut str = [0u8; 16];
        let mut stack: Vec<Frame> = Vec::with_capacity(64);
        for a in 0..4u8 {
            stack.push(Frame {
                p: self.set_intv(a),
                d: 1,
                c: a,
            });
        }
        while let Some(top) = stack.pop() {
            if top.d > 0 {
                str[(depth - top.d) as usize] = top.c;
            }
            let mut ok = [Sai::default(); 4];
            self.extend(&top.p, &mut ok, true);
            for a in 0..4u8 {
                str[(depth - top.d - 1) as usize] = a;
                if top.d != depth - 1 {
                    stack.push(Frame {
                        p: ok[a as usize],
                        d: top.d + 1,
                        c: a,
                    });
                } else {
                    let mut x = 0usize;
                    for (i, &b) in str.iter().take(depth as usize).enumerate() {
                        x |= (b as usize) << (i * 2);
                    }
                    s[x] = ok[a as usize];
                }
            }
        }
    }
}

impl Default for Bwt {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 朴素参照：text + '$' 的完整后缀数组（仅测试用）。
    fn naive_sa(text: &[u8]) -> Vec<u64> {
        let n = text.len();
        let mut sa: Vec<usize> = (0..=n).collect();
        sa.sort_by(|&a, &b| {
            // '$' 在末尾，视为最小
            for i in 0.. {
                let ca = if a + i == n { None } else { Some(text[a + i]) };
                let cb = if b + i == n { None } else { Some(text[b + i]) };
                match (ca, cb) {
                    (None, None) => return std::cmp::Ordering::Equal,
                    (None, Some(_)) => return std::cmp::Ordering::Less,
                    (Some(_), None) => return std::cmp::Ordering::Greater,
                    (Some(x), Some(y)) if x != y => return x.cmp(&y),
                    _ => {}
                }
            }
            std::cmp::Ordering::Equal
        });
        sa.iter().map(|&x| x as u64).collect()
    }

    /// 从 text（2-bit，不含 $）构建 Bwt（测试用，走 gen_sa 路径）。
    ///
    /// 完整 BWT：B[i] = T$[SA[i]-1]（环绕），其中 SA[i]==0 的行是 '$'；
    /// primary = S^{-1}(0)，即 '$' 所在行。$-removed BWT 去掉该行。
    fn build(text: &[u8], sa_bit: u32) -> Bwt {
        let sa = naive_sa(text);
        let n = text.len() as u64;
        let primary = sa.iter().position(|&x| x == 0).unwrap() as u64;
        let mut bwt_str = Vec::with_capacity(n as usize);
        for &s in &sa {
            if s == 0 {
                continue; // '$' 行（primary）
            }
            bwt_str.push(text[(s - 1) as usize]);
        }
        let mut bwt = Bwt::from_raw(&bwt_str, n, primary);
        bwt.gen_sa(sa_bit);
        bwt
    }

    #[test]
    fn rank_matches_naive() {
        let text: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 2, 2, 3, 3, 0, 0, 1, 2, 3];
        let bwt = build(&text, 2);
        // 重建原始文本验证 SA
        let n = text.len() as u64;
        let mut recovered = vec![0u8; n as usize];
        let mut j = bwt.primary;
        for i in 0..n {
            j = bwt.inv_psi(j);
            recovered[i as usize] = bwt.b0(j - u64::from(j > bwt.primary));
        }
        // Ψ 游走从 primary 出发按文本逆序重建
        assert_eq!(recovered, text.iter().copied().rev().collect::<Vec<u8>>());
        // rank11 对暴力。注意 k 是"含 $ 坐标"：$-removed 下标 i 对应 k = i + (i >= primary)
        let s: Vec<u8> = bwt_to_string(&bwt);
        for c in 0..4u8 {
            let mut cnt = 0u64;
            for i in 0..=n {
                let k = i + u64::from(i >= bwt.primary);
                assert_eq!(bwt.rank11(k, c), cnt, "rank11 i={i} k={k} c={c}");
                if i < n && s[i as usize] == c {
                    cnt += 1;
                }
            }
        }
        // rank1a 对暴力
        for i in 0..=n {
            let k = i + u64::from(i >= bwt.primary);
            let mut cnt = [0u64; 4];
            bwt.rank1a(k, &mut cnt);
            for c in 0..4u8 {
                assert_eq!(cnt[c as usize], bwt.rank11(k, c), "rank1a i={i} k={k} c={c}");
            }
        }
    }

    fn bwt_to_string(bwt: &Bwt) -> Vec<u8> {
        (0..bwt.seq_len).map(|k| bwt.b0(k)).collect()
    }

    #[test]
    fn sa_matches_naive() {
        let text: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 2, 2, 3, 3, 0, 0, 1, 2, 3];
        let bwt = build(&text, 2);
        let sa_full = naive_sa(&text);
        // sa(k) 的 k 为"含 $ 坐标"（0..=n）：直接对应完整 SA 的行号。
        // k=0 在 C 版中故意置 -1（区间坐标恒 ≥1，不会查询 0 行）；
        // k=primary 是 $ 所在行，任何 SA 区间都不会包含它。两者均跳过。
        assert_eq!(bwt.sa(0), u64::MAX);
        for k in 1..=text.len() as u64 {
            if k == bwt.primary {
                continue;
            }
            assert_eq!(bwt.sa(k), sa_full[k as usize], "with-$ row {k}");
        }
    }

    #[test]
    fn smem_finds_exact_match() {
        let text: Vec<u8> = vec![0, 1, 2, 3, 3, 2, 1, 0, 0, 1, 2, 3, 3, 2, 1, 0];
        let bwt = build(&text, 2);
        let q: Vec<u8> = vec![1, 2, 3, 3, 2]; // 出现在 [1,6)
        let mut p = Sai::default();
        let mut x = 0i64;
        let mut found = Vec::new();
        while x < q.len() as i64 {
            x = bwt.smem(&q, x, 3, 1, &mut p);
            if p.size > 0 {
                found.push((p.info >> 32, p.info & 0xffffffff, p.size));
            }
        }
        assert!(!found.is_empty());
        // 最长 SMEM 应覆盖整个查询
        assert!(found.iter().any(|(s, e, _)| *s == 0 && *e == 5));
    }

    #[test]
    fn dbg_bwt_string() {
        let text: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 2, 2, 3, 3, 0, 0, 1, 2, 3];
        let bwt = build(&text, 2);
        let sa = naive_sa(&text);
        let primary = sa.iter().position(|&x| x == 0).unwrap();
        // 朴素 $-removed BWT 串
        let expect: Vec<u8> = sa
            .iter()
            .enumerate()
            .filter(|&(row, _)| row != primary)
            .map(|(_, &s)| text[(s - 1) as usize])
            .collect();
        let got = bwt_to_string(&bwt);
        assert_eq!(got, expect, "BWT 串不一致");
        // rank11 直接对 BWT 串暴力
        let mut cnt = [0u64; 4];
        for i in 0..=bwt.seq_len {
            let k = i + u64::from(i >= bwt.primary);
            for c in 0..4u8 {
                assert_eq!(
                    bwt.rank11(k, c),
                    cnt[c as usize],
                    "rank11 i={i} k={k} c={c} bwt_str={got:?} primary={}",
                    bwt.primary
                );
            }
            if i < bwt.seq_len {
                cnt[got[i as usize] as usize] += 1;
            }
        }
    }

    #[test]
    fn save_load_roundtrip() {
        let text: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 2, 2, 3, 3, 0, 0];
        let bwt = build(&text, 2);
        let tmp = std::env::temp_dir().join("minibwa-rust-test.mbw");
        bwt.save(&tmp).unwrap();
        let bwt2 = Bwt::load(&tmp).unwrap();
        assert_eq!(bwt2.primary, bwt.primary);
        assert_eq!(bwt2.l2, bwt.l2);
        assert_eq!(bwt2.data, bwt.data);
        assert_eq!(bwt2.sa, bwt.sa);
        assert_eq!(bwt2.sa_bit, bwt.sa_bit);
    }

    #[test]
    fn kmer_cache_matches_extend() {
        let text: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 2, 2, 3, 3, 0, 0];
        let mut bwt = build(&text, 2);
        bwt.cache(3);
        let pre = bwt.pre.as_ref().unwrap();
        // 暴力：对每个 3-mer 用 extend 逐步搜索
        for (kmer, slot) in pre.iter().enumerate().take(64) {
            let bases = [(kmer >> 4) as u8 & 3, (kmer >> 2) as u8 & 3, kmer as u8 & 3];
            let mut p = bwt.set_intv(bases[0]);
            let mut ok = [Sai::default(); 4];
            for &c in &bases[1..] {
                bwt.extend(&p, &mut ok, true);
                p = ok[c as usize];
            }
            assert_eq!(*slot, p, "kmer {kmer:02b}");
        }
    }

    /// 对照 C 版 `mb_bwt_back`：cache 命中后游标必须停在 k-mer 左侧，
    /// SMEM 区间与无缓存路径一致。min_len-1 >= pre_len 才会走缓存。
    #[test]
    fn smem_with_cache_matches_without() {
        let text: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 2, 2, 3, 3, 0, 0];
        let bwt = build(&text, 2);
        let q = text.as_slice();
        let collect = |bwt: &Bwt| {
            let mut out = Vec::new();
            let mut p = Sai::default();
            let mut x = 0i64;
            while x < q.len() as i64 {
                x = bwt.smem(q, x, 5, 1, &mut p);
                if p.size > 0 {
                    out.push(p);
                }
            }
            out
        };
        let no_cache = collect(&bwt);
        let mut cached = bwt;
        cached.cache(3);
        assert_eq!(no_cache, collect(&cached));
    }
}
