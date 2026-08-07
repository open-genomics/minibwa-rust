// jemalloc：替换默认分配器，多线程场景显著提升
#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use anyhow::Result;
use clap::{Parser, Subcommand};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser, Debug)]
#[command(
    name = "minibwa-rust",
    author,
    version,
    about = "minibwa 的 Rust 重写：BWT 索引 + 变长种子 + minimap2 风格链化与比对",
    arg_required_else_help = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// 为参考基因组 FASTA 建索引（生成 <prefix>.l2b 与 <prefix>.mbw）
    Index {
        /// 参考基因组 FASTA（支持 .gz）
        fasta: String,
        /// 输出前缀（默认与输入文件相同）
        prefix: Option<String>,
        /// 模糊碱基随机种子
        #[arg(short = 's', long = "seed", default_value_t = 11u64)]
        seed: u64,
        /// SA 采样率 1/(1<<INT)
        #[arg(short = 'u', long = "sa-sample", default_value_t = 4u32)]
        sa_bit: u32,
        /// 线程数
        #[arg(short = 't', long = "threads", default_value_t = 4usize)]
        threads: usize,
    },
    /// 从 .l2b 还原参考基因组序列（调试用，与 C 版 getref 一致）
    Getref {
        /// .l2b 文件路径
        l2b: String,
    },
    /// 比对 reads（SAM/PAF 输出）
    Map {
        /// 索引前缀（<prefix>.l2b / <prefix>.mbw）
        reference: String,
        /// reads FASTA/FASTQ（支持 .gz）
        reads1: String,
        /// 第二端 reads（可选）
        reads2: Option<String>,
        /// 输出 PAF（默认 SAM）
        #[arg(short = 'f')]
        paf: bool,
        /// 只做链化不做碱基比对（隐含 -fP）
        #[arg(long = "chain-only")]
        chain_only: bool,
        /// 线程数
        #[arg(short = 't', default_value_t = 1)]
        threads: usize,
        /// 最小种子长度
        #[arg(short = 'k')]
        min_len: Option<i32>,
        /// 种子最大出现次数
        #[arg(short = 'c')]
        max_occ: Option<i32>,
        /// 最大 gap
        #[arg(short = 'g')]
        max_gap: Option<i32>,
        /// 带宽
        #[arg(short = 'w')]
        bw: Option<i32>,
        /// 长带宽
        #[arg(short = 'W')]
        bw_long: Option<i32>,
        /// 最小链得分
        #[arg(short = 'm')]
        min_chain_score: Option<i32>,
        /// 次优/最优得分比
        #[arg(short = 'p')]
        pri_ratio: Option<f32>,
        /// 最多保留次要比对数
        #[arg(short = 'N')]
        best_n: Option<i32>,
        /// adap 模式下短读长度上限
        #[arg(short = 'l')]
        max_sr_len: Option<i32>,
        /// 预设（sr/lr/adap）
        #[arg(short = 'x', default_value = "adap")]
        preset: String,
        /// 跳过配对
        #[arg(short = 'P')]
        no_pairing: bool,
        /// 取查询坐标最小的为 primary
        #[arg(short = '5')]
        primary5: bool,
        /// 不输出未比对 reads
        #[arg(short = 'u')]
        no_unmap: bool,
        /// 次要比对输出得分阈值（相对 best 的比例）
        #[arg(long = "outs")]
        out_s: Option<f32>,
        /// 最多输出的次要比对数
        #[arg(long = "outn")]
        out_n: Option<i32>,
        /// 输出文件（默认 stdout）
        #[arg(short = 'o')]
        out: Option<String>,
    },
    /// 兼容 bwa-mem 的旧式命令行（不推荐）
    #[command(disable_help_flag = true)]
    Mem {
        /// 索引前缀（<prefix>.l2b / <prefix>.mbw）
        reference: String,
        /// reads1（FASTA/FASTQ，支持 .gz）
        reads1: String,
        /// reads2（可选，配对端）
        reads2: Option<String>,
        /// 线程数
        #[arg(short = 't', default_value_t = 1)]
        threads: usize,
        /// 最小种子长度
        #[arg(short = 'k')]
        min_len: Option<i32>,
        /// 带宽
        #[arg(short = 'w')]
        bw: Option<i32>,
        /// 种子最大出现次数
        #[arg(short = 'c')]
        max_occ: Option<i32>,
        /// mask 等级
        #[arg(short = 'D')]
        mask_level: Option<f32>,
        /// 最小链得分
        #[arg(short = 'W')]
        min_chain_score: Option<i32>,
        /// mate rescue 候选数；0 跳过 rescue
        #[arg(short = 'm')]
        max_rescue: Option<i32>,
        /// 跳过 mate rescue
        #[arg(short = 'S')]
        no_rescue: bool,
        /// 跳过配对
        #[arg(short = 'P')]
        no_pairing: bool,
        /// 输出所有比对
        #[arg(short = 'a')]
        all: bool,
        /// max secondary alignments at the XA tag
        #[arg(short = 'h')]
        xa_max: Option<i32>,
        /// 跳过 XA 的得分阈值（相对 best 的比例）
        #[arg(short = 'z')]
        out_s: Option<f32>,
        /// SAM read group 行（如 '@RG\tID:foo\tSM:bar'）
        #[arg(short = 'R')]
        rg: Option<String>,
        /// 输出文件（默认 stdout）
        #[arg(short = 'o')]
        out: Option<String>,
        /// 取查询坐标最小的为 primary
        #[arg(short = '5')]
        primary5: bool,
        /// 抑制 DP 得分低于 INT×{-A} 的比对
        #[arg(short = 'T')]
        min_dp_max: Option<i32>,
        /// 复制 FASTA/Q 注释到输出
        #[arg(short = 'C')]
        copy_comment: bool,
        /// 补充比对使用软剪
        #[arg(short = 'Y')]
        supp_soft: bool,
        /// 不输出未比对 reads
        #[arg(short = 'u')]
        no_unmap: bool,
        /// insert size 分布 mean[,stddev[,max[,min]]]
        #[arg(short = 'I')]
        isize: Option<String>,
    },
    /// 调试：输出 SMEM 种子区间（与 C 版 fastmap 单序列模式一致）
    Fastmap {
        /// 索引前缀（<prefix>.l2b / <prefix>.mbw）
        prefix: String,
        /// reads FASTA/FASTQ（支持 .gz）
        reads: String,
        /// 最小种子长度
        #[arg(short = 'l', default_value_t = 19)]
        min_len: i32,
        /// 最小区间大小
        #[arg(short = 's', default_value_t = 1)]
        min_occ: i32,
        /// 输出坐标的最大区间大小
        #[arg(short = 'w', default_value_t = 20)]
        max_size_out: u64,
    },
    /// 打印版本号
    Version,
}

/// 一个 fragment（1 或 2 条 segment）。
#[derive(Debug, Clone, Copy)]
struct Fragment {
    off: usize, // 在 batch.seqs 中的起始下标
    cnt: usize, // 1 或 2
}

/// mini-batch：平坦的序列数组 + fragment 索引。
struct Batch {
    seqs: Vec<minibwa_rust::io::kseq::SeqRecord>,
    frags: Vec<Fragment>,
}

/// 去掉尾部 /[12] 后的 qname 长度（对应 `mb_qname_len`）。
fn qname_len(s: &str) -> usize {
    let b = s.as_bytes();
    let l = b.len();
    if l >= 3 && b[l - 1].is_ascii_digit() && b[l - 2] == b'/' {
        l - 2
    } else {
        l
    }
}

/// 两个 qname 是否为同一 fragment 的两端（对应 `mb_qname_same`）。
fn qname_same(s1: &str, s2: &str) -> bool {
    let l1 = qname_len(s1);
    let l2 = qname_len(s2);
    l1 == l2 && s1.as_bytes()[..l1] == s2.as_bytes()[..l2]
}

/// 单文件读取：按 qname 分组（相邻同 qname 的两条记录组成 fragment）。
fn read_batch(
    iter: &mut minibwa_rust::io::kseq::SeqIter<Box<dyn std::io::Read>>,
    mb_size: usize,
) -> Result<Batch> {
    let mut seqs = Vec::new();
    let mut frags = Vec::new();
    let mut bp = 0usize;
    let mut pending: Option<minibwa_rust::io::kseq::SeqRecord> = None;
    loop {
        if bp >= mb_size && pending.is_none() && !seqs.is_empty() {
            break;
        }
        let rec = match iter.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(e.into()),
            None => break,
        };
        bp += rec.seq.len();
        match pending.take() {
            None => pending = Some(rec),
            Some(p) => {
                if qname_same(&p.name, &rec.name) {
                    let off = seqs.len();
                    seqs.push(p);
                    seqs.push(rec);
                    frags.push(Fragment { off, cnt: 2 });
                } else {
                    let off = seqs.len();
                    seqs.push(p);
                    frags.push(Fragment { off, cnt: 1 });
                    pending = Some(rec);
                }
            }
        }
    }
    if let Some(p) = pending {
        let off = seqs.len();
        seqs.push(p);
        frags.push(Fragment { off, cnt: 1 });
    }
    Ok(Batch { seqs, frags })
}

/// 双文件读取：两个文件交替配对成 fragment。
fn read_batch_pe(
    iter1: &mut minibwa_rust::io::kseq::SeqIter<Box<dyn std::io::Read>>,
    iter2: &mut minibwa_rust::io::kseq::SeqIter<Box<dyn std::io::Read>>,
    mb_size: usize,
) -> Result<Batch> {
    let mut seqs = Vec::new();
    let mut frags = Vec::new();
    let mut bp = 0usize;
    loop {
        if bp >= mb_size && !seqs.is_empty() {
            break;
        }
        match (iter1.next(), iter2.next()) {
            (Some(Ok(a)), Some(Ok(b))) => {
                let off = seqs.len();
                bp += a.seq.len() + b.seq.len();
                seqs.push(a);
                seqs.push(b);
                frags.push(Fragment { off, cnt: 2 });
            }
            (None, None) => break,
            _ => {
                // 一端 EOF 另一端仍有记录
                eprintln!("[W] query files have different number of records; extra records skipped.");
                break;
            }
        }
    }
    Ok(Batch { seqs, frags })
}

/// 配对端 fragment 去掉尾部 "/[12]"（对应 C 版 step0 的 trim）。
fn trim_pair_suffix(seqs: &mut [minibwa_rust::io::kseq::SeqRecord], frags: &[Fragment]) {
    for f in frags {
        if f.cnt != 2 {
            continue;
        }
        let n0 = seqs[f.off].name.len();
        let n1 = seqs[f.off + 1].name.len();
        if n0 < 3 || n0 != n1 {
            continue;
        }
        let b0 = seqs[f.off].name.as_bytes();
        let b1 = seqs[f.off + 1].name.as_bytes();
        if b0[n0 - 1] != b1[n1 - 1] && b0[n0 - 2] == b'/' {
            seqs[f.off].name.truncate(n0 - 2);
            seqs[f.off + 1].name.truncate(n1 - 2);
        }
    }
}

/// 解析 -I mean[,stddev[,max[,min]]]（对应 C 版 `set_ins_size`）。
fn set_ins_size(opt: &mut minibwa_rust::opt::Opt, arg: &str) {
    let mut parts = arg.splitn(4, ',');
    let avg = parts.next().unwrap_or("").parse::<f64>().unwrap_or(0.0);
    let avg_i = (0.499 + avg) as i32;
    let std = match parts.next() {
        Some(s) if !s.is_empty() => (0.499 + s.parse::<f64>().unwrap_or(0.0)) as i32,
        _ => (avg_i as f64 * 0.1) as i32,
    };
    let (hi, lo) = match (parts.next(), parts.next()) {
        (Some(h), Some(l)) => (h.parse::<i64>().unwrap_or(0) as i32, l.parse::<i64>().unwrap_or(0) as i32),
        _ => (avg_i + std * 4, avg_i - std * 4),
    };
    opt.pe_avg = avg_i;
    opt.pe_std = std;
    opt.pe_hi = hi;
    opt.pe_lo = if lo < 1 { 1 } else { lo };
    opt.flag |= minibwa_rust::opt::flag::PE_PREDEF;
}

fn run_map(
    reference: &str,
    reads1: &str,
    reads2: Option<&str>,
    opt: &minibwa_rust::opt::Opt,
    out_path: Option<&str>,
    rg: Option<&str>,
) -> Result<()> {
    use minibwa_rust::format::{fmt_paf, fmt_sam, fmt_sam_hdr};
    use minibwa_rust::index::Idx;
    use minibwa_rust::io::kseq::SeqIter;
    use minibwa_rust::l2bit::open_fasta;
    use minibwa_rust::mapcore::map_one;
    use minibwa_rust::opt::flag;
    use minibwa_rust::pe::{pair, pes_predefined, pestat, Pestat};
    use rayon::prelude::*;
    use std::io::Write;

    let idx = Idx::load(reference)?;
    let mut out: Box<dyn Write> = match out_path {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };

    let is_paf = opt.flag & flag::PAF != 0;
    let mut rg_id = None;
    if !is_paf {
        let (hdr, rg) = fmt_sam_hdr(&idx.l2b, env!("CARGO_PKG_VERSION"), rg);
        rg_id = rg;
        out.write_all(hdr.as_bytes())?;
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opt.n_thread.max(1) as usize)
        .build()
        .map_err(|e| anyhow::anyhow!("初始化线程池失败: {e}"))?;
    let mb_size = opt.mb_size.max(1) as usize;

    let reader1 = open_fasta(reads1)?;
    let mut it1 = SeqIter::new(reader1);
    let mut it2 = reads2.map(|p| SeqIter::new(open_fasta(p).unwrap()));

    loop {
        let mut batch = if let Some(it2) = it2.as_mut() {
            read_batch_pe(&mut it1, it2, mb_size)?
        } else {
            read_batch(&mut it1, mb_size)?
        };
        if batch.seqs.is_empty() {
            break;
        }
        trim_pair_suffix(&mut batch.seqs, &batch.frags);

        // 1. 单端映射（rayon 并行）
        let mapped: Vec<minibwa_rust::mapcore::MapOut> = pool.install(|| {
            batch
                .seqs
                .par_iter()
                .map(|rec| map_one(opt, &idx, &rec.name, &rec.seq))
                .collect()
        });
        let mut seg_hits: Vec<Vec<minibwa_rust::hits::Hit>> =
            mapped.into_iter().map(|m| m.hits).collect();
        let mut n_hit: Vec<i32> = seg_hits.iter().map(|h| h.len() as i32).collect();

        // 2. 配对编排（pestat → pair）
        let n_frag = batch.frags.len();
        let n_pe = batch.frags.iter().filter(|f| f.cnt == 2).count();
        let is_pe = opt.flag & flag::PE != 0
            && n_frag < batch.seqs.len()
            && opt.flag & (flag::NO_PAIRING | flag::NO_ALN) == 0;
        if is_pe {
            let pes: [Pestat; 4] = if opt.flag & flag::PE_PREDEF != 0 || n_pe < 20 {
                pes_predefined(opt)
            } else {
                let seg_off: Vec<usize> = batch.frags.iter().map(|f| f.off).collect();
                let seg_cnt: Vec<usize> = batch.frags.iter().map(|f| f.cnt).collect();
                let mut p = [Pestat::default(); 4];
                pestat(opt, n_frag, &seg_off, &seg_cnt, &seg_hits, &mut p);
                p
            };
            // 逐 fragment 配对
            for f in &batch.frags {
                if f.cnt != 2 {
                    continue;
                }
                let off = f.off;
                let mut nh = [n_hit[off], n_hit[off + 1]];
                let mut hits2: [Vec<minibwa_rust::hits::Hit>; 2] = [
                    std::mem::take(&mut seg_hits[off]),
                    std::mem::take(&mut seg_hits[off + 1]),
                ];
                let qlen2 = [
                    batch.seqs[off].seq.len() as i32,
                    batch.seqs[off + 1].seq.len() as i32,
                ];
                let qseq2 = [batch.seqs[off].seq.clone(), batch.seqs[off + 1].seq.clone()];
                pair(opt, &idx.l2b, &mut nh, &mut hits2, &pes, &qlen2, &qseq2);
                n_hit[off] = nh[0];
                n_hit[off + 1] = nh[1];
                let [h0, h1] = hits2;
                seg_hits[off] = h0;
                seg_hits[off + 1] = h1;
            }
        }

        // 3. 按 fragment 输出
        let mut s = String::new();
        for f in &batch.frags {
            let seg_st = f.off;
            let seg_en = f.off + f.cnt;
            for i in seg_st..seg_en {
                let t = &batch.seqs[i];
                let mate_qlen = if seg_en - seg_st > 1 {
                    let mate_idx = if i != seg_en - 1 { i + 1 } else { seg_st };
                    batch.seqs[mate_idx].seq.len()
                } else {
                    0
                };
                let n_seg = seg_en - seg_st;
                let seg_idx = i - seg_st;
                let mate_pri = if n_seg > 1 {
                    let midx = if i != seg_en - 1 { i + 1 } else { seg_st };
                    minibwa_rust::format::get_sam_pri(&seg_hits[midx])
                } else {
                    None
                };
                let qual: Option<&[u8]> = if t.qual.is_empty() { None } else { Some(&t.qual) };
                let hits = &seg_hits[i];
                if !hits.is_empty() {
                    let mut n_sec = 0i32;
                    for (hi, h) in hits.iter().enumerate() {
                        let is_pri = h.parent == h.id;
                        if !is_pri && n_sec >= opt.out_n {
                            continue;
                        }
                        if !is_pri {
                            let p = &hits[h.parent as usize];
                            let p_score = p.p.as_ref().map_or(p.score as f64, |pp| pp.dp_max as f64);
                            let h_score = h.p.as_ref().map_or(h.score as f64, |pp| pp.dp_max as f64);
                            if h_score < opt.out_s as f64 * p_score {
                                continue;
                            }
                        }
                        if is_paf {
                            fmt_paf(&mut s, &idx.l2b, &t.name, t.seq.len(), Some(h), n_seg, seg_idx);
                        } else {
                            fmt_sam(
                                &mut s,
                                &idx.l2b,
                                &t.name,
                                &t.seq,
                                qual,
                                hits,
                                Some(hi),
                                n_seg,
                                seg_idx,
                                mate_pri,
                                mate_qlen,
                                rg_id.as_deref(),
                                opt,
                            );
                        }
                        if !is_pri {
                            n_sec += 1;
                        }
                    }
                } else if opt.flag & flag::NO_UNMAP == 0 {
                    if is_paf {
                        fmt_paf(&mut s, &idx.l2b, &t.name, t.seq.len(), None, n_seg, seg_idx);
                    } else {
                        fmt_sam(
                            &mut s,
                            &idx.l2b,
                            &t.name,
                            &t.seq,
                            qual,
                            hits,
                            None,
                            n_seg,
                            seg_idx,
                            mate_pri,
                            mate_qlen,
                            rg_id.as_deref(),
                            opt,
                        );
                    }
                }
            }
        }
        out.write_all(s.as_bytes())?;
    }
    out.flush()?;
    Ok(())
}

fn run_fastmap(
    prefix: &str,
    reads: &str,
    min_len: i32,
    min_occ: i32,
    max_size_out: u64,
) -> Result<()> {
    use minibwa_rust::bwt::Sai;
    use minibwa_rust::dna::NT4;
    use minibwa_rust::index::Idx;
    use minibwa_rust::io::kseq::SeqIter;
    use minibwa_rust::l2bit::open_fasta;
    use std::io::Write;

    let idx = Idx::load(prefix)?;
    let reader = open_fasta(reads)?;
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());

    for rec in SeqIter::new(reader) {
        let rec = rec?;
        let seq: Vec<u8> = rec.seq.iter().map(|&b| NT4[b as usize]).collect();
        writeln!(out, "SQ\t{}\t{}", rec.name, seq.len())?;
        // 与 C 版一致：fastmap 不建 kmer 缓存，直接走标准 SMEM
        let mut p = Sai::default();
        let mut x = 0i64;
        let mut intervals = Vec::new();
        while {
            x = idx.bwt.smem(&seq, x, min_len as i64, min_occ as i64, &mut p);
            if p.size > 0 {
                intervals.push(p);
            }
            x < seq.len() as i64
        } {}
        for p in &intervals {
            let st = p.info >> 32;
            let en = p.info & 0xffff_ffff;
            write!(out, "EM\t{st}\t{en}\t{}", p.size)?;
            let len = en - st;
            if p.size <= max_size_out {
                for row in p.x[0]..p.x[0] + p.size {
                    let sa = idx.bwt.sa(row);
                    let mut rev = false;
                    let (tid, cst) = idx.l2b.intv2cid(sa, sa + len, &mut rev);
                    if tid < 0 {
                        write!(out, "\t.")?;
                    } else {
                        let strand = if rev { '-' } else { '+' };
                        write!(out, "\t{}:{}{}", idx.l2b.ctg[tid as usize].name, strand, cst + 1)?;
                    }
                }
            } else {
                write!(out, "\t*")?;
            }
            writeln!(out)?;
        }
        writeln!(out, "//")?;
    }
    out.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Index {
            fasta,
            prefix,
            seed,
            sa_bit,
            threads,
        } => {
            let prefix = prefix.unwrap_or_else(|| fasta.clone());
            minibwa_rust::index::build(&fasta, &prefix, seed, sa_bit, threads)
        }
        Commands::Getref { l2b } => minibwa_rust::l2bit::getref(&l2b),
        Commands::Fastmap {
            prefix,
            reads,
            min_len,
            min_occ,
            max_size_out,
        } => run_fastmap(&prefix, &reads, min_len, min_occ, max_size_out),
        Commands::Map {
            reference,
            reads1,
            reads2,
            paf,
            chain_only,
            threads,
            min_len,
            max_occ,
            max_gap,
            bw,
            bw_long,
            min_chain_score,
            pri_ratio,
            best_n,
            max_sr_len,
            preset,
            no_pairing,
            primary5,
            no_unmap,
            out_s,
            out_n,
            out,
        } => {
            let mut opt = minibwa_rust::opt::Opt::default();
            if opt.apply_preset(&preset).is_err() {
                eprintln!("ERROR: 未知预设 '{preset}'");
                std::process::exit(1);
            }
            // 仅显式传入的参数才覆盖预设（与 C 版两遍解析语义一致）
            if let Some(v) = min_len {
                opt.min_len = v;
            }
            if let Some(v) = max_occ {
                opt.max_occ = v;
            }
            if let Some(v) = max_gap {
                opt.max_gap = v;
            }
            if let Some(v) = bw {
                opt.bw = v;
            }
            if let Some(v) = bw_long {
                opt.bw_long = v;
            }
            if let Some(v) = min_chain_score {
                opt.min_chain_score = v;
            }
            if let Some(v) = pri_ratio {
                opt.pri_ratio = v;
            }
            if let Some(v) = best_n {
                opt.best_n = v;
            }
            if let Some(v) = max_sr_len {
                opt.max_sr_len = v;
            }
            if let Some(v) = out_s {
                opt.out_s = v;
            }
            if let Some(v) = out_n {
                opt.out_n = v;
            }
            opt.n_thread = threads.max(1) as i32;
            if no_pairing {
                opt.flag |= minibwa_rust::opt::flag::NO_PAIRING;
            }
            if primary5 {
                opt.flag |= minibwa_rust::opt::flag::PRIMARY5;
            }
            if no_unmap {
                opt.flag |= minibwa_rust::opt::flag::NO_UNMAP;
            }
            if paf {
                opt.flag |= minibwa_rust::opt::flag::PAF;
            }
            if chain_only {
                opt.flag |= minibwa_rust::opt::flag::NO_ALN;
            }
            if opt.flag & minibwa_rust::opt::flag::NO_ALN != 0 {
                opt.flag |= minibwa_rust::opt::flag::NO_PAIRING | minibwa_rust::opt::flag::PAF;
            }
            run_map(&reference, &reads1, reads2.as_deref(), &opt, out.as_deref(), None)
        }
        Commands::Mem {
            reference,
            reads1,
            reads2,
            threads,
            min_len,
            bw,
            max_occ,
            mask_level,
            min_chain_score,
            max_rescue,
            no_rescue,
            no_pairing,
            all,
            xa_max,
            out_s,
            rg,
            out,
            primary5,
            min_dp_max,
            copy_comment,
            supp_soft,
            no_unmap,
            isize,
        } => {
            let mut opt = minibwa_rust::opt::Opt::default();
            opt.flag |= minibwa_rust::opt::flag::WRITE_MD; // bwa-mem 总是写 MD
            if let Some(v) = min_len {
                opt.min_len = v;
            }
            if let Some(v) = bw {
                opt.bw = v;
            }
            if let Some(v) = max_occ {
                opt.max_occ = v;
            }
            if let Some(v) = mask_level {
                opt.mask_level = v;
            }
            if let Some(v) = min_chain_score {
                opt.min_chain_score = v;
            }
            if let Some(v) = max_rescue {
                opt.max_rescue = v;
            }
            if no_rescue {
                opt.max_rescue = 0;
            }
            if no_pairing {
                opt.flag |= minibwa_rust::opt::flag::NO_PAIRING;
            }
            if all {
                opt.out_n = 1_000_000;
            }
            if let Some(v) = xa_max {
                opt.xa_max = v;
            }
            if let Some(v) = out_s {
                opt.out_s = v;
            }
            if primary5 {
                opt.flag |= minibwa_rust::opt::flag::PRIMARY5;
            }
            if let Some(v) = min_dp_max {
                opt.min_dp_max = v;
            }
            if copy_comment {
                opt.flag |= minibwa_rust::opt::flag::COPY_COMMENT;
            }
            if supp_soft {
                opt.flag |= minibwa_rust::opt::flag::SUPP_SOFT;
            }
            if no_unmap {
                opt.flag |= minibwa_rust::opt::flag::NO_UNMAP;
            }
            if let Some(v) = isize {
                set_ins_size(&mut opt, &v);
            }
            opt.n_thread = threads.max(1) as i32;
            run_map(&reference, &reads1, reads2.as_deref(), &opt, out.as_deref(), rg.as_deref())
        }
        Commands::Version => {
            println!("{VERSION}");
            Ok(())
        }
    }
}
