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
    /// 打印版本号
    Version,
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
