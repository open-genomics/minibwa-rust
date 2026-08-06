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
    Mem {
        reference: String,
        reads1: String,
        reads2: Option<String>,
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

fn run_map(
    reference: &str,
    reads1: &str,
    reads2: Option<&str>,
    opt: &minibwa_rust::opt::Opt,
    out_path: Option<&str>,
) -> Result<()> {
    use minibwa_rust::format::fmt_paf;
    use minibwa_rust::index::Idx;
    use minibwa_rust::io::kseq::SeqIter;
    use minibwa_rust::l2bit::open_fasta;
    use minibwa_rust::mapcore::map_one;
    use rayon::prelude::*;
    use std::io::Write;

    if reads2.is_some() {
        anyhow::bail!("配对端比对将在 M6 里程碑交付");
    }
    if opt.flag & minibwa_rust::opt::flag::PAF == 0 {
        anyhow::bail!("SAM 输出将在 M5 里程碑交付（当前可用 --chain-only/-f 输出 PAF）");
    }
    let idx = Idx::load(reference)?;

    let reader = open_fasta(reads1)?;
    let mut out: Box<dyn Write> = match out_path {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };

    // 按 mb_size 碱基为一批，rayon 并行处理（保持输入顺序输出）
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opt.n_thread.max(1) as usize)
        .build()
        .map_err(|e| anyhow::anyhow!("初始化线程池失败: {e}"))?;
    let mb_size = opt.mb_size.max(1) as usize;
    let mut batch: Vec<(String, Vec<u8>)> = Vec::new();
    let mut batch_bp = 0usize;
    let mut records = SeqIter::new(reader);
    loop {
        let mut eof = false;
        while batch_bp < mb_size {
            match records.next() {
                Some(Ok(rec)) => {
                    batch_bp += rec.seq.len();
                    batch.push((rec.name, rec.seq));
                }
                Some(Err(e)) => return Err(e.into()),
                None => {
                    eof = true;
                    break;
                }
            }
        }
        if batch.is_empty() {
            break;
        }
        let outputs: Vec<String> = pool.install(|| {
            batch
            .par_iter()
            .map(|(name, seq)| {
                let res = map_one(opt, &idx, name, seq);
                let mut s = String::new();
                if res.hits.is_empty() {
                    if opt.flag & minibwa_rust::opt::flag::NO_UNMAP == 0 {
                        fmt_paf(&mut s, &idx.l2b, name, seq.len(), None, 1, 0);
                    }
                } else {
                    let mut n_sec = 0i32;
                    for h in &res.hits {
                        let is_pri = h.parent == h.id;
                        if !is_pri && n_sec >= opt.out_n {
                            continue;
                        }
                        if !is_pri {
                            let p = &res.hits[h.parent as usize];
                            if (h.score as f64) < opt.out_s as f64 * p.score as f64 {
                                continue;
                            }
                        }
                        fmt_paf(&mut s, &idx.l2b, name, seq.len(), Some(h), 1, 0);
                        if !is_pri {
                            n_sec += 1;
                        }
                    }
                }
                s
            })
            .collect()
        });
        for s in outputs {
            out.write_all(s.as_bytes())?;
        }
        batch.clear();
        batch_bp = 0;
        if eof {
            break;
        }
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
            run_map(&reference, &reads1, reads2.as_deref(), &opt, out.as_deref())
        }
        Commands::Mem { .. } => {
            eprintln!("ERROR: mem 子命令（bwa-mem 兼容 CLI）尚未实现");
            std::process::exit(1);
        }
        Commands::Version => {
            println!("{VERSION}");
            Ok(())
        }
    }
}
