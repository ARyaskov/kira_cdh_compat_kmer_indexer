//! CD-HIT-compatible k-mer indexing (“CD-HIT-NG”) in modern Rust (edition 2024).
//!
//! New in 0.2:
//! - O(1) rolling reverse-complement for canonicalization
//! - Optional postings compression (delta+varint)
//! - Optional PLR sampling for sub-log locate
//! - Fast ASCII→2-bit mapping (LUT) + optional AVX2 path behind `simd_x86`
//! - Property tests via `proptest`
//!
//! Public API remains stable. The `Range<u64>` returned by [`KmerIndex::locate`]
//! is opaque and must be passed unmodified to [`KmerIndex::postings`].
//!
//! See README for on-disk format v2.

mod builder;
pub mod encode;
mod index;
mod io;
mod plr;
mod radix; // lightweight PGM-like sampling

pub use crate::index::LocateStrategy;
pub use crate::io::Compression;
pub use builder::{BuildConfig, build_kmer_index_async, build_kmer_index_sync};
pub use encode::{canonical, encode_kmer, revcomp};
pub use index::{KmerIndex, Posting};
pub use io::KmerIndexWriter;

pub use kira_cdh_compat_fastq_reader::ReaderOptions;
pub use std::ops::Range;
