use clap::Parser;
use kira_cdh_compat_kmer_indexer::Compression;
use kira_cdh_compat_kmer_indexer::LocateStrategy;
use kira_cdh_compat_kmer_indexer::*;
use std::path::PathBuf;

/// Build a `.kdx` k-mer index from FASTQ.
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input FASTQ(.gz) path
    #[arg(short, long)]
    input: PathBuf,

    /// Output `.kdx` path
    #[arg(short, long)]
    output: PathBuf,

    /// K-mer length (<= 32)
    #[arg(short = 'k', long)]
    k: usize,

    /// Bucket bits (B), default 12
    #[arg(short = 'B', long, default_value_t = 12)]
    bucket_bits: u8,

    /// Disable canonicalization
    #[arg(long, default_value_t = false)]
    no_canonical: bool,

    /// Include heads[] (forced true when compression != none)
    #[arg(long, default_value_t = false)]
    heads: bool,

    /// Minimum read length
    #[arg(long, default_value_t = 0)]
    min_read_len: usize,

    /// Threads (sync builder)
    #[arg(long)]
    threads: Option<usize>,

    /// Compression: none|dv (delta+varint)
    #[arg(long, default_value = "none")]
    compression: String,

    /// Locate strategy: bin|plr
    #[arg(long, default_value = "bin")]
    locate: String,

    /// PLR stride (valid for locate=plr)
    #[arg(long, default_value_t = 256)]
    plr_stride: usize,
}

fn parse_compression(s: &str) -> Compression {
    match s {
        "none" => Compression::None,
        "dv" | "delta" | "varint" => Compression::DeltaVarint,
        _ => Compression::None,
    }
}

fn parse_locate(s: &str) -> LocateStrategy {
    match s {
        "bin" | "binary" => LocateStrategy::BinarySearch,
        "plr" | "pgm" => LocateStrategy::PlrSampling,
        _ => LocateStrategy::BinarySearch,
    }
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let compression = parse_compression(&args.compression);
    let locate = parse_locate(&args.locate);

    let cfg = BuildConfig::default()
        .with_bucket_bits(args.bucket_bits)
        .canonical(!args.no_canonical)
        .with_heads(args.heads)
        .min_read_len(args.min_read_len)
        .compression(compression)
        .locate_strategy(locate)
        .plr_stride(args.plr_stride);

    let cfg = match args.threads {
        Some(n) => cfg.threads(n),
        None => cfg,
    };
    let idx = build_kmer_index_sync(&args.input, Default::default(), args.k, cfg)?;

    KmerIndexWriter::new(&idx).write_to(&args.output)?;
    eprintln!(
        "Built kdx: k={}, B={}, canonical={}, heads={}, compression={:?}, locate={:?}",
        idx.k(),
        args.bucket_bits,
        idx.is_canonical(),
        args.heads,
        compression,
        locate
    );

    Ok(())
}
