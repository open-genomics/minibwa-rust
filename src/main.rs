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
        /// 索引前缀（<prefix>.l2b / <prefix>.mbw）或直接给参考 FASTA
        reference: String,
        /// reads FASTA/FASTQ（支持 .gz）
        reads1: String,
        /// 第二端 reads（可选）
        reads2: Option<String>,
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
        Commands::Map { .. } | Commands::Mem { .. } => {
            eprintln!("ERROR: map 子命令尚未实现（将在 M5/M6 里程碑交付）");
            std::process::exit(1);
        }
        Commands::Version => {
            println!("{VERSION}");
            Ok(())
        }
    }
}
