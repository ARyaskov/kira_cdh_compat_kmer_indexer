# kira_cdh_compat_kmer_indexer

CD-HIT-compatible k-mer indexing (“CD-HIT-NG”) for FASTQ in modern Rust (edition 2024).

- **Alphabet:** A=00, C=01, G=10, T/U=11 (2 bits/base)
- **Canonicalization:** O(1) rolling reverse-complement (`min(fwd, rc)`)
- **Ambiguity:** windows containing IUPAC/`N` are skipped (rolling state resets)
- **Sharding:** top-B bits on **MSB-aligned** codes (1K–16K buckets typical)
- **Parallel build:** per-thread local buffers → deterministic count/prefix-sum/move
- **Sort:** in-bucket LSD radix sort (8 passes × 8 bits) on `u64` codes
- **On-disk:** single `.kdx` file, mmap-friendly; lazy decode for compressed postings
- **Locate:** binary search or PLR sampling (two-level) to narrow the search range
- **Memory target:** ~16 B/entry (8B code + 8B posting), plus optional heads/samples

> **CD-HIT compatibility:** encoding, canonicalization and `N`‐handling match CD-HIT semantics so the index can drive a CD-HIT-style pipeline.
> The `.kdx` file format is this crate’s own (compact & mmap-oriented).

---

## Quick Start

```bash
# Build CLI
cargo build --release

# Create a k-mer index (31-mers, 4096 buckets, canonicalized)
./target/release/kdx-build \
  --input reads.fastq.gz \
  --output reads.kdx \
  --k 31 \
  --bucket-bits 12
````

Reopen and query in your own code:

```rust
use kira_cdh_compat_kmer_indexing as kdx;
use std::path::Path;

let idx = kdx::KmerIndex::open_mmap(Path::new("reads.kdx"))?;
let code = {
    // 31-mer MSB-aligned, canonicalized
    let lsb = kdx::encode::encode_kmer(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
    kdx::encode::canonical(lsb, 31)
};
let range = idx.locate(code);
for p in idx.postings(range) {
    // p.read_id (u32), p.pos (u32)
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

---

## Installation

```toml
# Cargo.toml
[dependencies]
kira_cdh_compat_kmer_indexer = "*"
```

**Feature flags**

* `mmap` *(default)* — enables `memmap2` + `bytemuck` for mmap reads
* `async` — enables the async builder (requires `tokio`)
* `simd_x86` — enables optional AVX2 path for ASCII→2-bit mapping (auto-fallback)

---

## Library API

### Build (sync)

```rust
use kira_cdh_compat_kmer_indexer::*;
use kira_cdh_compat_kmer_indexer::io::Compression;
use kira_cdh_compat_kmer_indexer::index::LocateStrategy;
use std::path::Path;

let idx = build_kmer_index_sync(
    Path::new("reads.fastq.gz"),
    ReaderOptions::default(),
    31,
    BuildConfig::default()
        .with_bucket_bits(12)          // B ∈ [10..14] typical
        .canonical(true)               // default true
        .with_heads(true)              // required when compression != None
        .compression(Compression::DeltaVarint)
        .locate_strategy(LocateStrategy::PlrSampling)
        .plr_stride(256)
)?;

// Optional: serialize to disk (v2 format)
KmerIndexWriter::new(&idx).write_to(Path::new("reads.kdx"))?;
```

### Build (async)

Enable the `async` feature and use:

```rust
let idx = build_kmer_index_async(
    Path::new("reads.fastq.gz"),
    ReaderOptions::default(),
    31,
    BuildConfig::default()
).await?;
```

### Query

```rust
let kdx = KmerIndex::open_mmap(Path::new("reads.kdx"))?;
let msb_code = encode::canonical(encode::encode_kmer(b"ACGT...window").unwrap(), 31);
let range = kdx.locate(msb_code);     // opaque; pass it unchanged
let posts = kdx.postings(range);      // &[Posting { read_id, pos }]
```

> `range` is **opaque** (internally tagged with bucket id). Always pass it unchanged to `postings()`.

---

## CLI

```
kdx-build --input <FASTQ[.gz]> --output <reads.kdx> --k <1..=32>
          [--bucket-bits B] [--no-canonical] [--heads]
          [--compression {none|dv}] [--locate {bin|plr}] [--plr-stride N]
          [--min-read-len N] [--threads N]
```

**Recommended presets**

* Typical large datasets: `--k 31 --bucket-bits 12 --locate plr --plr-stride 256`
* Disk-sensitive output: `--compression dv --heads` (heads are required by `dv`)

---

## Design & Performance Notes

### MSB-aligned codes

* Internally we store 2k bits in the **upper** bits of `u64` (MSB-aligned).
  This makes bucket extraction trivial and keeps radix passes cache-friendly.

### Rolling canonicalization (O(1))

Per base:

```
fwd = ((fwd << 2) | v) & mask
rc  = (rc >> 2) | ((v ^ 0b11) << (2*(k-1)))
code = min(fwd, rc) << (64 - 2k)
```

Ambiguous base (IUPAC/`N`) → reset `fwd`, `rc`, and window length.

### Sharding

* Bucket id: **top-B** bits of the MSB code.
* We clamp `B` to `min(14, 2k)` to avoid degenerate bucketing for small `k`.

### Sorting

* LSD radix sort, 8×8-bit passes, stable via counting + prefix sums.
* Postings are permuted alongside codes.

### Locate strategies

* **BinarySearch** (default): classic lower/upper bound in `codes[bucket]`.
* **PlrSampling** (optional): store every *stride*-th `(key, pos)`; at query time, upper-bound in the samples and binary-search only within a narrow window. Typical `stride` 256–512.

### Postings compression (optional)

`Compression::DeltaVarint`:

* Requires `heads[]` (first index of each unique code run).
* For each run, postings are sorted by `(read_id, pos)` and encoded as LEB128: first pair absolute, subsequent as deltas.
* On mmap open, a bucket is **decompressed lazily** on first `postings()` access and cached in memory.

### Memory/throughput tips

* Build in `--release` with `RUSTFLAGS="-C target-cpu=native"`.
* Use `--threads` to cap parallelism or let Rayon use logical cores.
* `B=12` (4096 buckets) is a solid default for k=31.

---

## File Format (`.kdx`)

Version **2** (mmap-friendly, supports compression/samples). All integers are **little-endian**; payload aligned to 8 bytes.

### Header (v2)

| Field             | Type | Notes                         |
| ----------------- | ---- | ----------------------------- |
| `magic`           | u32  | `"KDX1"`                      |
| `version`         | u32  | `2`                           |
| `k`               | u16  | k-mer length                  |
| `canonical`       | u8   | 0/1                           |
| `bucket_bits`     | u8   | B                             |
| `buckets`         | u32  | `1 << B`                      |
| `total_entries`   | u64  | sum of all postings           |
| `dir_offset`      | u64  | byte offset of the directory  |
| `compression`     | u8   | 0=None, 1=DeltaVarint         |
| `locate_strategy` | u8   | 0=BinarySearch, 1=PlrSampling |
| `reserved0`       | u16  | 0                             |
| `reserved1`       | u32  | 0                             |

### Directory entry (per bucket, v2)

| Field         | Type | Meaning                                                                |
| ------------- | ---- | ---------------------------------------------------------------------- |
| `off_codes`   | u64  | byte offset of `codes[]` (u64 MSB-aligned keys)                        |
| `len_codes`   | u64  | number of `u64` codes                                                  |
| `off_posts`   | u64  | if compressed: **byte** offset of encoded postings; else: array offset |
| `len_posts`   | u64  | if compressed: **byte** length; else: number of `Posting` entries      |
| `total_posts` | u64  | total postings count (used when compressed)                            |
| `off_heads`   | u64  | byte offset of `heads[]` (`u32`), or 0                                 |
| `len_heads`   | u64  | number of `u32` head entries, or 0                                     |
| `off_samples` | u64  | byte offset of PLR sample pairs `(u64 key, u64 pos)`, or 0             |
| `len_samples` | u64  | number of sample pairs, or 0                                           |

> The crate can still **open** v1 files (no compression/samples) but **writes** v2.

---

## Determinism & Reproducibility

* Fixed read → per-thread shards → deterministic global merge.
* Stable radix passes; per-run posting sort is deterministic.
* Identical inputs and config produce byte-identical `.kdx`.

---

## Testing & Benchmarks

* Unit tests: `cargo test`
* Property tests (`proptest`): included for rolling encode consistency
* Micro-benchmarks (Criterion): `cargo bench` (enable in your workspace as needed)

---

## Configuration Cheatsheet

| Option                 | Default | Notes                                                          |
| ---------------------- | ------- | -------------------------------------------------------------- |
| `k`                    | —       | ≤ 32 (u64). `k>32` would need u128 or two u64 (future option). |
| `bucket_bits`          | 12      | Clamped to `min(14, 2k)`                                       |
| `canonical`            | true    | Rolling O(1) implementation                                    |
| `heads`                | false   | **Required** if `compression=DeltaVarint`                      |
| `compression`          | None    | `DeltaVarint` reduces on-disk size; lazy per-bucket decode     |
| `locate_strategy`      | Binary  | `PlrSampling` reduces locate time                              |
| `plr_stride`           | 256     | 256–512 is a good starting range                               |
| `min_read_len`         | 0       | Skip too-short reads early                                     |
| `threads` (sync build) | auto    | Override to pin build parallelism                              |

---

## FAQ

**Q:** Why MSB-aligned codes?
**A:** Faster bucket extraction and cache-friendly radix passes. The alignment is internal; API accepts/returns MSB-aligned codes for lookups.

**Q:** Is the index stable across runs?
**A:** Yes. Build is deterministic for identical inputs and config.

**Q:** Do I need `heads[]`?
**A:** Not for uncompressed postings. For `DeltaVarint` compression, `heads[]` is required to delineate code runs during (de)compression.

**Q:** How does PLR sampling compare to a full PGM index?
**A:** PLR sampling is a lightweight two-level approach: it gets you close to PGM locate speeds with minimal build time and footprint. A full PGM index will be added later.

---

## Versioning

* `0.2.x` — rolling RC (O(1)), optional `DeltaVarint` compression, PLR sampling, `.kdx` v2.
* `0.1.x` — baseline implementation, `.kdx` v1.

---

## License

GPLv2

