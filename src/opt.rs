//! 比对参数，对应 C 版 minibwa 的 `mb_opt_t`（minibwa.h / options.c）。

/// 输出与行为标志位（对应 `MB_F_*`）。
pub mod flag {
    pub const PAF: u64 = 0x1; // 以 PAF 格式输出
    pub const NO_UNMAP: u64 = 0x2; // 不输出未比对序列
    pub const COPY_COMMENT: u64 = 0x4; // FASTX comment 复制到输出
    pub const PE: u64 = 0x8; // 配对端模式
    pub const LONG: u64 = 0x10; // 长序列模式
    pub const EQX: u64 = 0x20; // CIGAR 使用 =/X
    pub const NO_KALLOC: u64 = 0x40; // Rust 版无 kalloc，恒置位
    pub const NO_ALN: u64 = 0x80; // 跳过碱基级比对
    pub const PE_PREDEF: u64 = 0x100; // 使用预定义的 PE 分布
    pub const WRITE_DS: u64 = 0x200; // 输出 ds:Z
    pub const WRITE_CS: u64 = 0x400; // 输出 cs:Z
    pub const WRITE_MD: u64 = 0x800; // 输出 MD:Z
    pub const SECOND_SEQ: u64 = 0x1000; // SAM 中为次要比对输出 SEQ
    pub const SUPP_SOFT: u64 = 0x2000; // 补充比对使用软剪
    pub const ADAP: u64 = 0x4000; // 自适应模式
    pub const PRIMARY5: u64 = 0x8000; // Hi-C 用
    pub const NO_PAIRING: u64 = 0x10000; // 不配对
    pub const METH: u64 = 0x20000; // 甲基化模式
}

/// 比对参数（字段命名与 C 版保持一致以便对照）。
#[derive(Debug, Clone)]
pub struct Opt {
    pub flag: u64,
    // 种子
    pub min_len: i32,     // 最小种子长度
    pub max_sub_occ: i32, // SMEM 出现次数低于该值时寻找更短种子
    pub max_occ: i32,     // 区间最大出现次数
    // 通用算法
    pub bw: i32, // 带宽
    pub bw_long: i32,
    pub max_gap: i32,    // 链中断的最大空窗
    pub max_sr_len: i32, // adap 模式下超过该长度按长读处理
    // 链化
    pub max_chain_skip: i32,
    pub max_chain_iter: i32,
    pub min_chain_score: i32,
    pub chain_gap_scale: f32,
    // hit 处理
    pub mask_level: f32,
    pub mask_len: i32,
    pub pri_ratio: f32,
    pub best_n: i32,
    // 比对打分
    pub a: i32,      // 匹配
    pub b: i32,      // 错配
    pub b_ts: i32,   // 转换错配
    pub b_ambi: i32, // 模糊碱基错配
    pub q: i32,      // 开 gap
    pub q2: i32,     // 开长 gap
    pub e: i32,      // 延伸 gap
    pub e2: i32,     // 延伸长 gap
    pub end_bonus: i32,
    pub min_dp_max: i32, // 最低得分阈值（×a）
    pub zdrop: i32,
    pub zdrop_inv: i32,
    pub min_ksw_len: i32,
    // 配对
    pub max_pe_ins: i32,
    pub max_rescue: i32,
    pub pen_unpair: i32,
    pub pe_avg: i32,
    pub pe_std: i32,
    pub pe_lo: i32,
    pub pe_hi: i32,
    // I/O
    pub sb_len: i32, // 批量 SMEM 的碱基总数
    pub sb_seq: i32, // 批量 SMEM 的序列条数
    pub n_thread: i32,
    pub out_n: i32, // 最多输出的次要比对数
    pub out_s: f32,
    pub seed: i32,
    pub xa_max: i32,
    pub mb_size: i64, // mini-batch 大小
    pub max_mb_size: i64,
    pub max_sw_mat: i64,
    pub cap_kalloc: i64,
}

impl Default for Opt {
    /// 对应 `mb_opt_init()`：先 reset 再应用 adap 预设。
    fn default() -> Self {
        let mut opt = Self::reset();
        opt.apply_preset("adap").expect("adap preset 恒有效");
        opt
    }
}

impl Opt {
    /// 对应 `mb_opt_reset()`：基础默认值（不含预设）。
    pub fn reset() -> Self {
        Opt {
            flag: flag::NO_KALLOC,
            min_len: 19,
            max_sub_occ: 10,
            max_occ: 250,
            bw: 0,
            bw_long: 20000,
            max_gap: 0,
            max_sr_len: 325,
            max_chain_skip: 25,
            max_chain_iter: 5000,
            min_chain_score: 0,
            chain_gap_scale: 0.8,
            mask_level: 0.5,
            mask_len: i32::MAX,
            pri_ratio: 0.5,
            best_n: 0,
            a: 2,
            b: 8,
            b_ts: 0,
            b_ambi: 1,
            q: 12,
            q2: 23,
            e: 2,
            e2: 1,
            end_bonus: 0,
            min_dp_max: 0,
            zdrop: 0,
            zdrop_inv: 0,
            min_ksw_len: 0,
            max_pe_ins: 10000,
            max_rescue: 10,
            pen_unpair: 17,
            pe_avg: 400,
            pe_std: 100,
            pe_lo: 50,
            pe_hi: 800,
            sb_len: 1_000_000,
            sb_seq: 24,
            n_thread: 1,
            out_n: 0,
            out_s: 0.8,
            seed: 11,
            xa_max: 5,
            mb_size: 0,
            max_mb_size: 1_000_000_000,
            max_sw_mat: 100_000_000,
            cap_kalloc: 1 << 28,
        }
    }

    /// 对应 `mb_opt_preset()`。preset ∈ {"sr", "lr", "adap"}。
    pub fn apply_preset(&mut self, preset: &str) -> Result<(), String> {
        *self = Self::reset();
        match preset {
            "sr" | "adap" => {
                self.flag |= flag::PE;
                if preset == "adap" {
                    self.flag |= flag::ADAP;
                }
                self.min_dp_max = 30;
                self.flag |= flag::ADAP; // C 版在 sr 分支也置 ADAP（保持原样）
                self.bw = 100;
                self.max_gap = 100;
                self.zdrop = 80;
                self.zdrop_inv = 80;
                self.best_n = 50;
                self.end_bonus = 10;
                self.min_chain_score = 25;
                self.min_ksw_len = 20;
                self.mb_size = 100_000_000;
            }
            "lr" => {
                self.flag |= flag::LONG;
                self.flag &= !flag::PE;
                self.min_dp_max = 50;
                self.bw = 500;
                self.max_gap = 5000;
                self.zdrop = 400;
                self.zdrop_inv = 240;
                self.best_n = 5;
                self.end_bonus = -1;
                self.min_chain_score = 40;
                self.min_ksw_len = 200;
                self.mb_size = 500_000_000;
            }
            _ => return Err(format!("unknown preset: {preset}")),
        }
        Ok(())
    }

    /// 是否为短读模式（对应 `mb_is_sr_mode`）。
    pub fn is_sr_mode(&self, qlen: i32) -> bool {
        !(self.flag & flag::LONG != 0 || (self.flag & flag::ADAP != 0 && qlen > self.max_sr_len))
    }

    /// 对应 `mb_opt_adap()`：按读长自适应调整参数（原地修改副本）。
    pub fn adap(&self, len: i32) -> Opt {
        const MIN_LEN: f64 = 100.0;
        const MID_LEN: f64 = 2000.0;
        let mut opt = self.clone();
        if self.flag & flag::ADAP == 0 {
            return opt;
        }
        let a = -(0.5f64.ln()) / (MID_LEN - MIN_LEN);
        let l = if f64::from(len) > MIN_LEN {
            f64::from(len)
        } else {
            MIN_LEN
        };
        let b = (-a * (l - MIN_LEN)).exp();
        if self.max_gap < 5000 {
            opt.max_gap = (5000.0 - f64::from(5000 - self.max_gap) * b + 0.499) as i32;
            if opt.max_gap > len {
                opt.max_gap = len;
            }
        }
        if self.bw < 500 {
            opt.bw = (500.0 - f64::from(500 - self.bw) * b + 0.499) as i32;
        }
        if opt.bw_long > len * 5 {
            opt.bw_long = len * 5;
        }
        if opt.bw_long < opt.bw {
            opt.bw_long = opt.bw;
        }
        if self.zdrop < 400 {
            opt.zdrop = (400.0 - f64::from(400 - self.zdrop) * b + 0.499) as i32;
        }
        if self.zdrop_inv < 240 {
            opt.zdrop_inv = (240.0 - f64::from(240 - self.zdrop_inv) * b + 0.499) as i32;
        }
        if self.best_n > 5 {
            opt.best_n = (f64::from(self.best_n - 5) * b + 5.0 + 0.499) as i32;
        }
        if self.min_dp_max < 50 {
            opt.min_dp_max = (50.0 - f64::from(50 - self.min_dp_max) * b + 0.499) as i32;
        }
        if self.min_chain_score < 40 {
            opt.min_chain_score = (40.0 - f64::from(40 - self.min_chain_score) * b + 0.499) as i32;
        }
        opt
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_adap_preset() {
        let opt = Opt::default();
        assert_eq!(opt.min_len, 19);
        assert_eq!(opt.bw, 100);
        assert_eq!(opt.zdrop, 80);
        assert_eq!(opt.a, 2);
        assert_eq!(opt.b, 8);
        assert!(opt.flag & flag::PE != 0);
        assert!(opt.flag & flag::ADAP != 0);
    }

    #[test]
    fn adap_interpolates_with_read_length() {
        let opt = Opt::default();
        // 短读：参数不变
        let o100 = opt.adap(100);
        assert_eq!(o100.bw, 100);
        // 长读：向 lr 端点靠拢
        let o2000 = opt.adap(2000);
        assert_eq!(o2000.bw, 300); // 500 - (500-100)*0.5 = 300
        let o_long = opt.adap(10000);
        assert!(o_long.bw > 450);
    }
}
