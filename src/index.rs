//! KmerIndex: in-memory (owned) and mmap-backed index with per-bucket layout.

use bytemuck::{Pod, Zeroable};
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, OnceLock};

use crate::encode::bucket_id_from_msb;
use crate::io::{BucketDir, Compression, FileHeader, KDX_MAGIC, KDX_VERSION2};
use crate::plr::PlrSamples;
use crate::radix::radix_sort_pairs_u64;
use thiserror::Error;

const RANGE_TAG_SHIFT: u32 = 48; // store bucket id in upper bits of start/end

/// Posting entry: `(read_id, pos)`.
#[repr(C)]
#[derive(Copy, Clone, Default, Pod, Zeroable, PartialEq, Eq, Debug)]
pub struct Posting {
    /// Read ordinal id (0-based).
    pub read_id: u32,
    /// Position within the read (0-based).
    pub pos: u32,
}

#[derive(Debug, Error)]
/// Errors returned by KmerIndex.
pub enum IndexError {
    /// I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// Invalid file format.
    #[error("Invalid KDX file: {0}")]
    Format(String),
    /// Bytemuck cast failed.
    #[error("Cast error: {0}")]
    Cast(String),
}

/// Locate acceleration strategy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocateStrategy {
    /// Plain binary search within a bucket.
    BinarySearch,
    /// Two-level sampling (PLR-like).
    PlrSampling,
}

/// Owned buckets data.
struct OwnedBuckets {
    codes: Vec<Vec<u64>>,
    posts: Vec<Vec<Posting>>,
    heads: Option<Vec<Vec<u32>>>,
    // Optional samples per bucket
    samples: Option<Vec<Option<PlrSamples>>>,
    strategy: LocateStrategy,
}

/// Mmap-backed data.
struct MmapData {
    map: memmap2::Mmap, // keep mmap alive
    header: FileHeader,
    dir: Vec<BucketDir>,
    // Lazy decompression cache for postings per bucket (only if compressed)
    post_cache: Vec<OnceLock<Arc<Vec<Posting>>>>,
    // Sampling tables per bucket (loaded from file when present)
    samples_cache: Vec<OnceLock<Option<PlrSamples>>>,
}

/// KmerIndex can be either owned or mmap-backed.
enum Storage {
    Owned(OwnedBuckets),
    Mmap(MmapData),
}

/// K-mer index: query-time API for both owned and mmap-backed layouts.
pub struct KmerIndex {
    k: usize,
    canonical: bool,
    bucket_bits: u8,
    buckets: u32,
    total_entries: u64,
    compression: Compression,
    storage: Storage,
}

impl KmerIndex {
    /// Create an empty owned index with configured geometry.
    pub(crate) fn new_owned(
        k: usize,
        canonical: bool,
        bucket_bits: u8,
        buckets: u32,
        with_heads: bool,
        compression: Compression,
        strategy: LocateStrategy,
    ) -> Self {
        let heads = if with_heads {
            Some(vec![Vec::<u32>::new(); buckets as usize])
        } else {
            None
        };
        KmerIndex {
            k,
            canonical,
            bucket_bits,
            buckets,
            total_entries: 0,
            compression,
            storage: Storage::Owned(OwnedBuckets {
                codes: vec![Vec::<u64>::new(); buckets as usize],
                posts: vec![Vec::<Posting>::new(); buckets as usize],
                heads,
                samples: match strategy {
                    LocateStrategy::PlrSampling => Some(vec![None; buckets as usize]),
                    _ => None,
                },
                strategy,
            }),
        }
    }

    /// Open a `.kdx` file via mmap (v1 or v2).
    pub fn open_mmap(path: &Path) -> Result<Self, IndexError> {
        use std::fs::File;
        use std::io::{Seek, SeekFrom};

        let file = File::open(path)?;
        let mut reader = std::io::BufReader::new(file);
        let header = FileHeader::read_from(&mut reader)?;

        if header.magic != KDX_MAGIC {
            return Err(IndexError::Format("bad magic".into()));
        }
        if header.version != KDX_VERSION2 && header.version != 1 {
            return Err(IndexError::Format("unsupported version".into()));
        }

        // Read directory
        reader.seek(SeekFrom::Start(header.dir_offset))?;
        let mut dir = Vec::with_capacity(header.buckets as usize);
        for _ in 0..header.buckets {
            dir.push(BucketDir::read_from(&mut reader, header.version)?);
        }

        // Map file
        let file = reader.into_inner();
        let map = unsafe { memmap2::MmapOptions::new().map(&file)? };

        let compression = if header.version == 1 {
            Compression::None
        } else {
            header.compression()
        };
        let strategy = if header.version == 1 {
            LocateStrategy::BinarySearch
        } else if header.locate_strategy_u8 == 1 {
            LocateStrategy::PlrSampling
        } else {
            LocateStrategy::BinarySearch
        };

        let post_cache = (0..header.buckets).map(|_| OnceLock::new()).collect();
        let samples_cache = (0..header.buckets).map(|_| OnceLock::new()).collect();

        Ok(KmerIndex {
            k: header.k as usize,
            canonical: header.canonical != 0,
            bucket_bits: header.bucket_bits,
            buckets: header.buckets,
            total_entries: header.total_entries,
            compression,
            storage: Storage::Mmap(MmapData {
                map,
                header,
                dir,
                post_cache,
                samples_cache,
            }),
        })
    }

    /// Bucket id by top-B bits of an MSB-aligned code.
    #[inline]
    pub fn bucket_id(&self, code_msb: u64) -> usize {
        bucket_id_from_msb(code_msb, self.bucket_bits)
    }

    /// Return the k used by this index.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Whether canonicalization was used.
    #[inline]
    pub fn is_canonical(&self) -> bool {
        self.canonical
    }

    /// Locate the postings range for a given MSB-aligned code.
    ///
    /// The returned range is **opaque** and must be passed as-is to [`postings`].
    pub fn locate(&self, code_msb: u64) -> Range<u64> {
        let b = self.bucket_id(code_msb);
        match &self.storage {
            Storage::Owned(o) => {
                let codes = &o.codes[b];
                if codes.is_empty() {
                    return tag_range(b as u32, 0..0);
                }
                let (start, end) = match o.strategy {
                    LocateStrategy::BinarySearch => lower_upper_bound(codes, code_msb),
                    LocateStrategy::PlrSampling => {
                        if let Some(Some(s)) = o.samples.as_ref().map(|v| &v[b]) {
                            let (lo, hi) = s.range_hint(code_msb, codes.len());
                            let (s2, e2) = lower_upper_bound_in(&codes[lo..hi], code_msb);
                            (lo + s2, lo + e2)
                        } else {
                            lower_upper_bound(codes, code_msb)
                        }
                    }
                };
                tag_range(b as u32, (start as u64)..(end as u64))
            }
            Storage::Mmap(m) => {
                let d = &m.dir[b];
                let codes = d.codes_slice::<u64>(&m.map).expect("codes slice");
                if codes.is_empty() {
                    return tag_range(b as u32, 0..0);
                }

                let (start, end) =
                    if (m.header.version >= KDX_VERSION2) && (m.header.locate_strategy_u8 == 1) {
                        let s_opt = m.samples_cache[b]
                            .get_or_init(|| d.samples_slice(&m.map))
                            .clone();
                        if let Some(s) = s_opt {
                            let (lo, hi) = s.range_hint(code_msb, codes.len());
                            let (s2, e2) = lower_upper_bound_in(&codes[lo..hi], code_msb);
                            (lo + s2, lo + e2)
                        } else {
                            lower_upper_bound(codes, code_msb)
                        }
                    } else {
                        lower_upper_bound(codes, code_msb)
                    };
                tag_range(b as u32, (start as u64)..(end as u64))
            }
        }
    }

    /// Borrow postings slice for the opaque `range` returned by [`locate`].
    ///
    /// If the on-disk file is compressed, the corresponding bucket is decompressed
    /// on first access and cached.
    pub fn postings(&self, range: Range<u64>) -> &[Posting] {
        let (b, local) = untag_range(range);
        match &self.storage {
            Storage::Owned(o) => &o.posts[b][local.start as usize..local.end as usize],
            Storage::Mmap(m) => {
                let cache = m.post_cache[b].get_or_init(|| {
                    // materialize postings for this bucket
                    let d = &m.dir[b];
                    let posts: Arc<Vec<Posting>> = match m.header.compression() {
                        Compression::None => {
                            let slice = d.posts_slice::<Posting>(&m.map).expect("posts");
                            Arc::new(slice.to_vec())
                        }
                        Compression::DeltaVarint => {
                            // heads required:
                            let heads = d.heads_slice::<u32>(&m.map).expect("heads");
                            let bytes = d.posts_bytes(&m.map).expect("posts bytes");
                            let mut out = Vec::<Posting>::with_capacity(d.total_posts as usize);
                            crate::io::decode_delta_varint(bytes, heads, &mut out).expect("decode");
                            Arc::new(out)
                        }
                    };
                    posts
                });
                &cache[local.start as usize..local.end as usize]
            }
        }
    }

    // -------- Internal helpers used by builder / writer --------

    pub(crate) fn buckets(&self) -> u32 {
        self.buckets
    }
    pub(crate) fn bucket_bits(&self) -> u8 {
        self.bucket_bits
    }
    pub(crate) fn total_entries_mut_add(&mut self, delta: u64) {
        self.total_entries = self.total_entries.saturating_add(delta);
    }
    pub(crate) fn owned_mut(&mut self) -> Option<&mut OwnedBuckets> {
        match &mut self.storage {
            Storage::Owned(o) => Some(o),
            _ => None,
        }
    }

    pub(crate) fn finalize_bucket(
        codes: &mut Vec<u64>,
        posts: &mut Vec<Posting>,
        heads: Option<&mut Vec<u32>>,
        build_samples: Option<&mut Option<PlrSamples>>,
        stride: usize,
    ) {
        if codes.is_empty() {
            return;
        }
        radix_sort_pairs_u64(codes, posts);
        if let Some(h) = heads {
            h.clear();
            let mut prev = codes[0];
            h.push(0u32);
            for (i, &c) in codes.iter().enumerate().skip(1) {
                if c != prev {
                    h.push(i as u32);
                    prev = c;
                }
            }
        }
        if let Some(slot) = build_samples {
            *slot = PlrSamples::build(codes, stride);
        }
    }

    pub(crate) fn storage_ref(&self) -> StorageRef<'_> {
        match &self.storage {
            Storage::Owned(o) => StorageRef::Owned {
                codes: &o.codes,
                posts: &o.posts,
                heads: o.heads.as_ref().map(|h| &h[..]),
                samples: o.samples.as_ref().map(|v| &v[..]),
                strategy: o.strategy,
            },
            Storage::Mmap(_) => StorageRef::Mmap,
        }
    }

    pub(crate) fn bucket_mut(
        &mut self,
        b: usize,
    ) -> (
        &mut Vec<u64>,
        &mut Vec<Posting>,
        Option<&mut Vec<u32>>,
        Option<&mut Option<PlrSamples>>,
    ) {
        match &mut self.storage {
            Storage::Owned(o) => {
                let codes = &mut o.codes[b];
                let posts = &mut o.posts[b];
                let heads = o.heads.as_mut().map(|h| &mut h[b]);
                let samples = o.samples.as_mut().map(|s| &mut s[b]);
                (codes, posts, heads, samples)
            }
            _ => unreachable!("mutation only valid for owned storage"),
        }
    }

    pub(crate) fn total_entries(&self) -> u64 {
        self.total_entries
    }
    pub(crate) fn compression(&self) -> Compression {
        self.compression
    }
    pub(crate) fn strategy(&self) -> LocateStrategy {
        match &self.storage {
            Storage::Owned(o) => o.strategy,
            Storage::Mmap(m) => {
                if m.header.version >= KDX_VERSION2 && m.header.locate_strategy_u8 == 1 {
                    LocateStrategy::PlrSampling
                } else {
                    LocateStrategy::BinarySearch
                }
            }
        }
    }
}

// ----- Opaque range tagging -----

#[inline]
fn tag_range(b: u32, r: Range<u64>) -> Range<u64> {
    let tag = (b as u64) << RANGE_TAG_SHIFT;
    (r.start | tag)..(r.end | tag)
}

#[inline]
fn untag_range(r: Range<u64>) -> (usize, Range<u64>) {
    let b = (r.start >> RANGE_TAG_SHIFT) as usize;
    let mask = (1u64 << RANGE_TAG_SHIFT) - 1;
    let local = (r.start & mask)..(r.end & mask);
    (b, local)
}

// ----- lower/upper bound helpers -----

fn lower_upper_bound(codes: &[u64], key: u64) -> (usize, usize) {
    lower_upper_bound_in(codes, key)
}

fn lower_upper_bound_in(codes: &[u64], key: u64) -> (usize, usize) {
    // lower_bound
    let mut lo = 0usize;
    let mut hi = codes.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if codes[mid] < key {
            lo = mid + 1
        } else {
            hi = mid
        }
    }
    let start = lo;
    // upper_bound
    hi = codes.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if codes[mid] <= key {
            lo = mid + 1
        } else {
            hi = mid
        }
    }
    (start, lo)
}

// Internal reference for writer
pub(crate) enum StorageRef<'a> {
    Owned {
        codes: &'a [Vec<u64>],
        posts: &'a [Vec<Posting>],
        heads: Option<&'a [Vec<u32>]>,
        samples: Option<&'a [Option<PlrSamples>]>,
        strategy: LocateStrategy,
    },
    Mmap,
}
