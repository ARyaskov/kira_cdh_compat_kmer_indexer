//! Builders: synchronous (rayon) and asynchronous (tokio) pipelines.
//! Adds O(1) rolling RC, optional compression, and optional PLR samples.

use crate::ReaderOptions;
use crate::encode::{MAP_LUT, align_msb, bucket_id_from_msb};
use crate::index::{KmerIndex, LocateStrategy, Posting};
use crate::io::{Compression, KmerIndexWriter};
use crate::plr::PlrSamples;
use crate::radix::radix_sort_pairs_u64;

use rayon::prelude::*;
use std::path::Path;

/// Build-time configuration.
#[derive(Clone)]
pub struct BuildConfig {
    bucket_bits: u8,
    canonical: bool,
    with_heads: bool,
    min_read_len: usize,
    threads: Option<usize>,
    async_buf_size: usize,
    compression: Compression,
    locate_strategy: LocateStrategy,
    plr_stride: usize,
}

impl Default for BuildConfig {
    fn default() -> Self {
        Self {
            bucket_bits: 12,
            canonical: true,
            with_heads: false,
            min_read_len: 0,
            threads: None,
            async_buf_size: 64 * 1024,
            compression: Compression::None,
            locate_strategy: LocateStrategy::BinarySearch,
            plr_stride: 256,
        }
    }
}

impl BuildConfig {
    /// Set number of top bits `B` used for bucketing (recommended 10..14).
    pub fn with_bucket_bits(mut self, b: u8) -> Self {
        self.bucket_bits = b;
        self
    }
    /// Enable/disable canonicalization (default: true).
    pub fn canonical(mut self, yes: bool) -> Self {
        self.canonical = yes;
        self
    }
    /// Include heads[] arrays (default: false).
    pub fn with_heads(mut self, yes: bool) -> Self {
        self.with_heads = yes;
        self
    }
    /// Minimum read length to consider (shorter reads are skipped).
    pub fn min_read_len(mut self, n: usize) -> Self {
        self.min_read_len = n;
        self
    }
    /// Fix the number of threads used by rayon (sync build).
    pub fn threads(mut self, n: usize) -> Self {
        self.threads = Some(n);
        self
    }
    /// Async batch size (records per batch).
    pub fn async_buf_size(mut self, n: usize) -> Self {
        self.async_buf_size = n;
        self
    }
    /// Compression strategy for postings.
    pub fn compression(mut self, c: Compression) -> Self {
        self.compression = c;
        self
    }
    /// Locate strategy.
    pub fn locate_strategy(mut self, s: LocateStrategy) -> Self {
        self.locate_strategy = s;
        self
    }
    /// PLR stride (only used when locate_strategy=PlrSampling).
    pub fn plr_stride(mut self, s: usize) -> Self {
        self.plr_stride = s.max(32);
        self
    }

    pub(crate) fn bucket_bits(&self) -> u8 {
        self.bucket_bits
    }
    pub(crate) fn canonical_flag(&self) -> bool {
        self.canonical
    }
    /// Effective heads requirement (true if user wants heads OR compression needs it).
    pub(crate) fn heads_effective(&self) -> bool {
        self.with_heads || matches!(self.compression, Compression::DeltaVarint)
    }
}

/// Helper: map any error to `std::io::Error` with kind=Other.
fn io_other<E: std::fmt::Display>(e: E) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::Other, format!("{e}"))
}

/// Build a k-mer index synchronously using rayon.
pub fn build_kmer_index_sync(
    path: &Path,
    reader_opts: ReaderOptions,
    k: usize,
    cfg: BuildConfig,
) -> Result<KmerIndex, std::io::Error> {
    use kira_cdh_compat_fastq_reader::FastqReader;

    assert!(k > 0 && k <= 32, "k must be 1..=32");

    let mut B = cfg.bucket_bits().min(14);
    let max_b = (2 * k) as u8;
    if B > max_b {
        B = max_b;
        eprintln!("[kdx] Adjusted bucket_bits to {} (<= 2k)", B);
    }

    let buckets = 1u32 << B;
    let mut idx = KmerIndex::new_owned(
        k,
        cfg.canonical_flag(),
        B,
        buckets,
        cfg.heads_effective(),
        cfg.compression,
        cfg.locate_strategy,
    );

    if let Some(n) = cfg.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    let mut reader = FastqReader::from_path(path, reader_opts).map_err(io_other)?;

    let bcount = buckets as usize;
    let block_size = 16_384usize;
    let mut read_id: u64 = 0;
    let mut total_entries = 0u64;

    loop {
        let mut block: Vec<(u64, Vec<u8>)> = Vec::with_capacity(block_size);
        let mut took_any = false;

        for item in (&mut reader).take(block_size) {
            took_any = true;
            let rec = item.map_err(io_other)?;
     
            if rec.seq.len() < cfg.min_read_len {
                read_id += 1;
                continue;
            }
            block.push((read_id, rec.seq.to_vec()));
            read_id += 1;
        }
        if !took_any {
            break;
        }
        if block.is_empty() {
            continue;
        }

        let shards: Vec<(Vec<Vec<u64>>, Vec<Vec<Posting>>)> = block
            .par_iter()
            .map(|(rid, seq)| {
                let mut codes: Vec<Vec<u64>> = (0..bcount).map(|_| Vec::new()).collect();
                let mut posts: Vec<Vec<Posting>> = (0..bcount).map(|_| Vec::new()).collect();
                extract_kmers_rolling(
                    &mut codes,
                    &mut posts,
                    *rid as u32,
                    seq,
                    k,
                    B,
                    cfg.canonical,
                );
                (codes, posts)
            })
            .collect();

        for (codes_s, posts_s) in shards {
            for b in 0..bcount {
                let (codes_b, posts_b, _heads_b, _samples_b) = idx.bucket_mut(b);
                codes_b.extend_from_slice(&codes_s[b]);
                posts_b.extend_from_slice(&posts_s[b]);
                total_entries += codes_s[b].len() as u64;
            }
        }
    }

    // Finalize buckets: sort, heads, samples
    let stride = cfg.plr_stride;
    for b in 0..(buckets as usize) {
        let (codes_b, posts_b, heads_b, samples_b) = idx.bucket_mut(b);
        KmerIndex::finalize_bucket(codes_b, posts_b, heads_b, samples_b, stride);
    }
    idx.total_entries_mut_add(total_entries);

    Ok(idx)
}

#[cfg(feature = "async")]
pub async fn build_kmer_index_async(
    path: &Path,
    reader_opts: ReaderOptions,
    k: usize,
    cfg: BuildConfig,
) -> Result<KmerIndex, std::io::Error> {
    use kira_cdh_compat_fastq_reader::AsyncFastqReader;

    assert!(k > 0 && k <= 32, "k must be 1..=32");

    let mut B = cfg.bucket_bits().min(14);
    let max_b = (2 * k) as u8;
    if B > max_b {
        B = max_b;
        eprintln!("[kdx] Adjusted bucket_bits to {} (<= 2k)", B);
    }

    let buckets = 1u32 << B;
    let mut idx = KmerIndex::new_owned(
        k,
        cfg.canonical_flag(),
        B,
        buckets,
        cfg.heads_effective(),
        cfg.compression,
        cfg.locate_strategy,
    );

    let mut ar = AsyncFastqReader::from_path(path, reader_opts)
        .await
        .map_err(io_other)?;
    let mut read_id: u64 = 0;
    let bcount = buckets as usize;
    let mut total_entries = 0u64;

    loop {
        let mut batch: Vec<(u64, Vec<u8>)> = Vec::with_capacity(cfg.async_buf_size);
        for _ in 0..cfg.async_buf_size {
            match ar.next_record().await {
                Some(item) => {
                    let rec = item.map_err(io_other)?;
                    if rec.seq.len() < cfg.min_read_len {
                        read_id += 1;
                        continue;
                    }
                    batch.push((read_id, rec.seq.to_vec()));
                    read_id += 1;
                }
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }

        let shards: Vec<(Vec<Vec<u64>>, Vec<Vec<Posting>>)> =
            tokio::task::spawn_blocking(move || {
                batch
                    .par_iter()
                    .map(|(rid, seq)| {
                        let mut codes: Vec<Vec<u64>> = (0..bcount).map(|_| Vec::new()).collect();
                        let mut posts: Vec<Vec<Posting>> =
                            (0..bcount).map(|_| Vec::new()).collect();
                        extract_kmers_rolling(
                            &mut codes,
                            &mut posts,
                            *rid as u32,
                            seq,
                            k,
                            B,
                            cfg.canonical,
                        );
                        (codes, posts)
                    })
                    .collect()
            })
            .await
            .expect("spawn_blocking failed");

        for (codes_s, posts_s) in shards {
            for b in 0..bcount {
                let (codes_b, posts_b, _heads_b, _samples_b) = idx.bucket_mut(b);
                codes_b.extend_from_slice(&codes_s[b]);
                posts_b.extend_from_slice(&posts_s[b]);
                total_entries += codes_s[b].len() as u64;
            }
        }
    }

    let stride = cfg.plr_stride;
    for b in 0..(buckets as usize) {
        let (codes_b, posts_b, heads_b, samples_b) = idx.bucket_mut(b);
        KmerIndex::finalize_bucket(codes_b, posts_b, heads_b, samples_b, stride);
    }
    idx.total_entries_mut_add(total_entries);

    Ok(idx)
}

#[cfg(not(feature = "async"))]
pub async fn build_kmer_index_async(
    _path: &Path,
    _reader_opts: ReaderOptions,
    _k: usize,
    _cfg: BuildConfig,
) -> Result<KmerIndex, std::io::Error> {
    panic!("build_kmer_index_async requires feature \"async\"");
}

// ---- Rolling extraction (O(1) canonicalization) ----

fn extract_kmers_rolling(
    codes: &mut [Vec<u64>],
    posts: &mut [Vec<Posting>],
    read_id: u32,
    seq: &[u8],
    k: usize,
    bbits: u8,
    canonicalize: bool,
) {
    if seq.len() < k {
        return;
    }

    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut fwd: u64 = 0;
    let mut rc: u64 = 0;
    let mut len: usize = 0;

    for (i, &b) in seq.iter().enumerate() {
        let v_raw = MAP_LUT[b as usize];
        if v_raw > 3 {
            // ambiguous: reset
            fwd = 0;
            rc = 0;
            len = 0;
            continue;
        }
        let v = v_raw as u64;

        // Update rolling fwd and rc
        fwd = ((fwd << 2) | v) & mask;
        rc = (rc >> 2) | ((v ^ 0b11) << (2 * (k - 1)));
        len += 1;

        if len >= k {
            let start_pos = (i + 1 - k) as u32;
            let code_lsb = fwd;
            let code_msb = if canonicalize {
                let c = if fwd <= rc { fwd } else { rc };
                align_msb(c, k)
            } else {
                align_msb(code_lsb, k)
            };
            let b = bucket_id_from_msb(code_msb, bbits);
            codes[b].push(code_msb);
            posts[b].push(Posting {
                read_id,
                pos: start_pos,
            });
        }
    }
}
