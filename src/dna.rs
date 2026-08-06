//! DNA 编码表与工具函数。
//!
//! 编码表与 C 版 minibwa 的 `kommon.c` 完全一致：
//! A/C/G/T → 0/1/2/3，其余碱基 → 4。

/// ASCII → 2-bit 编码（对应 `kom_nt4_table`）。A=0, C=1, G=2, T=3, 其余=4。
pub const NT4: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

/// ASCII 互补表（对应 `kom_comp_table`）：仅对 ACGT/acgt 及 IUPAC 简并码有效，
/// 其余字节原样保留。
pub const COMP: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    t[0] = 3;
    t[1] = 2;
    t[2] = 1;
    t[3] = 0;
    t[b'A' as usize] = b'T';
    t[b'C' as usize] = b'G';
    t[b'G' as usize] = b'C';
    t[b'T' as usize] = b'A';
    t[b'a' as usize] = b't';
    t[b'c' as usize] = b'g';
    t[b'g' as usize] = b'c';
    t[b't' as usize] = b'a';
    // IUPAC 简并码（与 kom_comp_table 一致）
    t[b'R' as usize] = b'Y';
    t[b'Y' as usize] = b'R';
    t[b'S' as usize] = b'S';
    t[b'W' as usize] = b'W';
    t[b'K' as usize] = b'M';
    t[b'M' as usize] = b'K';
    t[b'B' as usize] = b'V';
    t[b'V' as usize] = b'B';
    t[b'D' as usize] = b'H';
    t[b'H' as usize] = b'D';
    t[b'r' as usize] = b'y';
    t[b'y' as usize] = b'r';
    t[b's' as usize] = b's';
    t[b'w' as usize] = b'w';
    t[b'k' as usize] = b'm';
    t[b'm' as usize] = b'k';
    t[b'b' as usize] = b'v';
    t[b'v' as usize] = b'b';
    t[b'd' as usize] = b'h';
    t[b'h' as usize] = b'd';
    t
};

/// splitmix64（对应 `kom_splitmix64`）。用于给模糊碱基随机分配低 2 bit。
pub fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

/// 原地反向互补 ASCII 序列（对应 `kom_revcomp`）。
pub fn revcomp(seq: &mut [u8]) {
    let len = seq.len();
    for i in 0..len / 2 {
        let t = seq[len - i - 1];
        seq[len - i - 1] = COMP[seq[i] as usize];
        seq[i] = COMP[t as usize];
    }
    if len & 1 == 1 {
        seq[len / 2] = COMP[seq[len / 2] as usize];
    }
}

/// 2-bit 编码序列的反向互补：反转并取 3-x（对应 `mb_seq_rev` + 3-c 语义）。
pub fn revcomp_2bit(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = 3 - *b;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nt4_matches_kom_nt4_table() {
        assert_eq!(NT4[b'A' as usize], 0);
        assert_eq!(NT4[b'C' as usize], 1);
        assert_eq!(NT4[b'G' as usize], 2);
        assert_eq!(NT4[b'T' as usize], 3);
        assert_eq!(NT4[b'N' as usize], 4);
        assert_eq!(NT4[b'n' as usize], 4);
        assert_eq!(NT4[b'-' as usize], 4);
    }

    #[test]
    fn splitmix64_first_values() {
        // 与 C 版 kom_splitmix64(seed=11) 的前两次输出对照
        let mut x = 11u64;
        let a = splitmix64(&mut x);
        let b = splitmix64(&mut x);
        assert_eq!(a, 0x50f5647d2380309d);
        assert_ne!(a, b);
    }

    #[test]
    fn revcomp_ascii() {
        let mut s = b"ACGTN".to_vec();
        revcomp(&mut s);
        assert_eq!(&s, b"NACGT");
    }
}
