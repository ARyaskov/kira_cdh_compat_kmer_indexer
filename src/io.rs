//! On-disk `.kdx` format v2: header + directory + payload.
//! Sections are 8-byte aligned. All integers are little-endian.
//!
//! v1 (backward-compat open): no compression, no samples, older dir layout.
//! v2: adds compression + samples offsets, locate strategy.

use bytemuck::{Pod, Zeroable};
use byteorder::{LittleEndian as LE, ReadBytesExt, WriteBytesExt};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;

use crate::index::{IndexError, KmerIndex, LocateStrategy, Posting};
use crate::plr::PlrSamples;

pub const KDX_MAGIC: u32 = 0x4B_44_58_31; // "KDX1"
pub const KDX_VERSION2: u32 = 2;

#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Compression {
    /// Uncompressed `Posting` array (8 bytes per entry).
    None = 0,
    /// Per-run delta+varint (LEB128), requires heads[] and per-run sort (read_id,pos).
    DeltaVarint = 1,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct FileHeader {
    pub magic: u32,
    pub version: u32,
    pub k: u16,
    pub canonical: u8,
    pub bucket_bits: u8,
    pub buckets: u32,
    pub total_entries: u64,
    pub dir_offset: u64,
    // v2 fields:
    pub compression_u8: u8,
    pub locate_strategy_u8: u8,
    pub reserved0: u16,
    pub reserved1: u32,
}

impl FileHeader {
    pub fn write_to<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        w.write_u32::<LE>(self.magic)?;
        w.write_u32::<LE>(self.version)?;
        w.write_u16::<LE>(self.k)?;
        w.write_u8(self.canonical)?;
        w.write_u8(self.bucket_bits)?;
        w.write_u32::<LE>(self.buckets)?;
        w.write_u64::<LE>(self.total_entries)?;
        w.write_u64::<LE>(self.dir_offset)?;
        w.write_u8(self.compression_u8)?;
        w.write_u8(self.locate_strategy_u8)?;
        w.write_u16::<LE>(self.reserved0)?;
        w.write_u32::<LE>(self.reserved1)?;
        Ok(())
    }

    pub fn read_from<R: Read>(r: &mut R) -> std::io::Result<Self> {
        let magic = r.read_u32::<LE>()?;
        let version = r.read_u32::<LE>()?;
        let k = r.read_u16::<LE>()?;
        let canonical = r.read_u8()?;
        let bucket_bits = r.read_u8()?;
        let buckets = r.read_u32::<LE>()?;
        let total_entries = r.read_u64::<LE>()?;
        let dir_offset = r.read_u64::<LE>()?;
        let (compression_u8, locate_strategy_u8, reserved0, reserved1) = if version >= KDX_VERSION2
        {
            (
                r.read_u8()?,
                r.read_u8()?,
                r.read_u16::<LE>()?,
                r.read_u32::<LE>()?,
            )
        } else {
            (0u8, 0u8, 0u16, 0u32)
        };

        Ok(FileHeader {
            magic,
            version,
            k,
            canonical,
            bucket_bits,
            buckets,
            total_entries,
            dir_offset,
            compression_u8,
            locate_strategy_u8,
            reserved0,
            reserved1,
        })
    }

    pub fn compression(&self) -> Compression {
        match self.compression_u8 {
            0 => Compression::None,
            1 => Compression::DeltaVarint,
            _ => Compression::None,
        }
    }
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct BucketDir {
    pub off_codes: u64,
    pub len_codes: u64, // number of u64 codes
    // posts region:
    pub off_posts: u64,   // if compressed: byte offset; else: array offset
    pub len_posts: u64,   // if compressed: byte length; else: number of Postings
    pub total_posts: u64, // total postings count (for compressed decode)
    // heads region:
    pub off_heads: u64,
    pub len_heads: u64, // number of u32 heads entries (or 0)
    // samples region (optional):
    pub off_samples: u64,
    pub len_samples: u64, // number of sample pairs
}

impl BucketDir {
    pub fn write_to<W: Write>(&self, w: &mut W, version: u32) -> std::io::Result<()> {
        w.write_u64::<LE>(self.off_codes)?;
        w.write_u64::<LE>(self.len_codes)?;
        w.write_u64::<LE>(self.off_posts)?;
        w.write_u64::<LE>(self.len_posts)?;
        if version >= KDX_VERSION2 {
            w.write_u64::<LE>(self.total_posts)?;
            w.write_u64::<LE>(self.off_heads)?;
            w.write_u64::<LE>(self.len_heads)?;
            w.write_u64::<LE>(self.off_samples)?;
            w.write_u64::<LE>(self.len_samples)?;
        } else {
            // v1 layout (no total_posts/off_samples/len_samples)
            w.write_u64::<LE>(self.off_heads)?;
            w.write_u64::<LE>(self.len_heads)?;
        }
        Ok(())
    }

    pub fn read_from<R: Read>(r: &mut R, version: u32) -> std::io::Result<Self> {
        let off_codes = r.read_u64::<LE>()?;
        let len_codes = r.read_u64::<LE>()?;
        let off_posts = r.read_u64::<LE>()?;
        let len_posts = r.read_u64::<LE>()?;
        if version >= KDX_VERSION2 {
            let total_posts = r.read_u64::<LE>()?;
            let off_heads = r.read_u64::<LE>()?;
            let len_heads = r.read_u64::<LE>()?;
            let off_samples = r.read_u64::<LE>()?;
            let len_samples = r.read_u64::<LE>()?;
            Ok(BucketDir {
                off_codes,
                len_codes,
                off_posts,
                len_posts,
                total_posts,
                off_heads,
                len_heads,
                off_samples,
                len_samples,
            })
        } else {
            let off_heads = r.read_u64::<LE>()?;
            let len_heads = r.read_u64::<LE>()?;
            Ok(BucketDir {
                off_codes,
                len_codes,
                off_posts,
                len_posts,
                total_posts: len_posts,
                off_heads,
                len_heads,
                off_samples: 0,
                len_samples: 0,
            })
        }
    }

    pub fn codes_slice<'a, T: bytemuck::Pod>(
        &self,
        map: &'a memmap2::Mmap,
    ) -> Result<&'a [T], IndexError> {
        if self.len_codes == 0 {
            return Ok(&[]);
        }
        let start = self.off_codes as usize;
        let end = start + self.len_codes as usize * std::mem::size_of::<T>();
        let bytes = &map[start..end];
        bytemuck::try_cast_slice(bytes).map_err(|e| IndexError::Cast(format!("{e:?}")))
    }

    pub fn posts_slice<'a, T: bytemuck::Pod>(
        &self,
        map: &'a memmap2::Mmap,
    ) -> Result<&'a [T], IndexError> {
        if self.len_posts == 0 {
            return Ok(&[]);
        }
        let start = self.off_posts as usize;
        let end = start + self.len_posts as usize * std::mem::size_of::<T>();
        let bytes = &map[start..end];
        bytemuck::try_cast_slice(bytes).map_err(|e| IndexError::Cast(format!("{e:?}")))
    }

    pub fn posts_bytes<'a>(&self, map: &'a memmap2::Mmap) -> Result<&'a [u8], IndexError> {
        if self.len_posts == 0 {
            return Ok(&[]);
        }
        let start = self.off_posts as usize;
        let end = start + self.len_posts as usize;
        Ok(&map[start..end])
    }

    pub fn heads_slice<'a, T: bytemuck::Pod>(
        &self,
        map: &'a memmap2::Mmap,
    ) -> Result<&'a [T], IndexError> {
        if self.len_heads == 0 {
            return Ok(&[]);
        }
        let start = self.off_heads as usize;
        let end = start + self.len_heads as usize * std::mem::size_of::<T>();
        let bytes = &map[start..end];
        bytemuck::try_cast_slice(bytes).map_err(|e| IndexError::Cast(format!("{e:?}")))
    }

    pub fn samples_slice(&self, map: &memmap2::Mmap) -> Option<PlrSamples> {
        if self.len_samples == 0 {
            return None;
        }
        let start = self.off_samples as usize;
        let end = start + self.len_samples as usize * std::mem::size_of::<SamplePair>();
        let bytes = &map[start..end];
        let raw: &[SamplePair] = bytemuck::try_cast_slice(bytes).ok()?;
        let mut keys = Vec::with_capacity(raw.len());
        let mut pos = Vec::with_capacity(raw.len());
        for sp in raw {
            keys.push(sp.key);
            pos.push(sp.pos);
        }
        // infer stride
        let stride = if pos.len() >= 2 {
            (pos[1] - pos[0]).max(1) as usize
        } else {
            256
        };
        Some(PlrSamples { stride, keys, pos })
    }
}

/// Writer that serializes an in-memory [`KmerIndex`] to a `.kdx` file.
pub struct KmerIndexWriter<'a> {
    idx: &'a KmerIndex,
}

impl<'a> KmerIndexWriter<'a> {
    /// Create a new writer from an in-memory index.
    pub fn new(idx: &'a KmerIndex) -> Self {
        Self { idx }
    }

    /// Serialize to disk. The produced file is memmap-friendly and deterministic.
    pub fn write_to(&self, path: &Path) -> Result<(), IndexError> {
        use std::fs::OpenOptions;

        let mut file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)?;

        // We only support writing from owned storage.
        let (codes, posts, heads_opt, samples_opt, strategy) = match self.idx.storage_ref() {
            crate::index::StorageRef::Owned {
                codes,
                posts,
                heads,
                samples,
                strategy,
            } => (codes, posts, heads, samples, strategy),
            crate::index::StorageRef::Mmap => {
                return Err(IndexError::Format(
                    "writing from mmap-backed index is not implemented".into(),
                ));
            }
        };

        let buckets = self.idx.buckets();
        let compression = self.idx.compression();
        let mut header = FileHeader {
            magic: KDX_MAGIC,
            version: KDX_VERSION2,
            k: self.idx.k() as u16,
            canonical: if self.idx.is_canonical() { 1 } else { 0 },
            bucket_bits: self.idx.bucket_bits(),
            buckets,
            total_entries: self.idx.total_entries(),
            dir_offset: 0,
            compression_u8: compression as u8,
            locate_strategy_u8: match strategy {
                LocateStrategy::BinarySearch => 0,
                LocateStrategy::PlrSampling => 1,
            },
            reserved0: 0,
            reserved1: 0,
        };

        // Reserve space: header + dir
        header.write_to(&mut file)?;
        let dir_off = file.stream_position()?;
        let mut dir_vec = vec![BucketDir::default(); buckets as usize];
        for d in &dir_vec {
            d.write_to(&mut file, KDX_VERSION2)?;
        }

        // Payload
        let align8 = |pos: u64| (pos + 7) & !7;
        let mut write_pos = file.stream_position()?;

        for b in 0..(buckets as usize) {
            let d = &mut dir_vec[b];

            // codes
            let cs = &codes[b];
            if !cs.is_empty() {
                write_pos = align8(write_pos);
                file.seek(SeekFrom::Start(write_pos))?;
                d.off_codes = write_pos;
                d.len_codes = cs.len() as u64;
                let bytes = bytemuck::cast_slice::<u64, u8>(cs);
                file.write_all(bytes)?;
                write_pos += (cs.len() * std::mem::size_of::<u64>()) as u64;
            }

            // heads (optional, or required if compressed)
            if let Some(heads) = heads_opt {
                let hs = &heads[b];
                if !hs.is_empty() {
                    write_pos = align8(write_pos);
                    file.seek(SeekFrom::Start(write_pos))?;
                    d.off_heads = write_pos;
                    d.len_heads = hs.len() as u64;
                    let bytes = bytemuck::cast_slice::<u32, u8>(hs);
                    file.write_all(bytes)?;
                    write_pos += (hs.len() * std::mem::size_of::<u32>()) as u64;
                }
            }

            // samples (optional)
            if let Some(samples) = samples_opt {
                if let Some(Some(s)) = samples.get(b) {
                    let mut tmp = Vec::<SamplePair>::with_capacity(s.keys.len());
                    for i in 0..s.keys.len() {
                        tmp.push(SamplePair {
                            key: s.keys[i],
                            pos: s.pos[i],
                        });
                    }
                    write_pos = align8(write_pos);
                    file.seek(SeekFrom::Start(write_pos))?;
                    d.off_samples = write_pos;
                    d.len_samples = tmp.len() as u64;
                    let bytes = bytemuck::cast_slice::<SamplePair, u8>(&tmp);
                    file.write_all(bytes)?;
                    write_pos += (tmp.len() * std::mem::size_of::<SamplePair>()) as u64;
                }
            }

            // postings
            let ps = &posts[b];
            match compression {
                Compression::None => {
                    if !ps.is_empty() {
                        write_pos = align8(write_pos);
                        file.seek(SeekFrom::Start(write_pos))?;
                        d.off_posts = write_pos;
                        d.len_posts = ps.len() as u64; // elements
                        d.total_posts = ps.len() as u64;
                        let bytes = bytemuck::cast_slice::<Posting, u8>(ps);
                        file.write_all(bytes)?;
                        write_pos += (ps.len() * std::mem::size_of::<Posting>()) as u64;
                    }
                }
                Compression::DeltaVarint => {
                    write_pos = align8(write_pos);
                    file.seek(SeekFrom::Start(write_pos))?;
                    d.off_posts = write_pos;
                    d.total_posts = ps.len() as u64;
                    // heads required:
                    let heads = heads_opt.expect("heads required for compressed postings");
                    let hs = &heads[b];
                    let bytes_written = encode_delta_varint(ps, hs, &mut file)?;
                    d.len_posts = bytes_written as u64; // bytes length
                    write_pos += bytes_written as u64;
                }
            }
        }

        // Write directory
        let _end_pos = align8(write_pos);
        file.seek(SeekFrom::Start(dir_off))?;
        for d in &dir_vec {
            d.write_to(&mut file, KDX_VERSION2)?;
        }

        // Patch header with dir offset
        let mut header2 = header;
        header2.dir_offset = dir_off;
        file.seek(SeekFrom::Start(0))?;
        header2.write_to(&mut file)?;

        Ok(())
    }
}

// ---------------- Varint (LEB128) utils ----------------

#[inline]
fn write_varu32<W: Write>(mut x: u32, w: &mut W) -> std::io::Result<usize> {
    let mut n = 0usize;
    while x >= 0x80 {
        w.write_all(&[((x as u8) | 0x80)])?;
        x >>= 7;
        n += 1;
    }
    w.write_all(&[x as u8])?;
    Ok(n + 1)
}

#[inline]
fn read_varu32<R: Read>(r: &mut R) -> std::io::Result<u32> {
    let mut x: u32 = 0;
    let mut s = 0u32;
    loop {
        let b = r.read_u8()?;
        x |= ((b & 0x7F) as u32) << s;
        if (b & 0x80) == 0 {
            break;
        }
        s += 7;
    }
    Ok(x)
}

/// Encode postings per-run as delta+varint.
/// Assumes `heads` marks the start index of each unique code run.
pub fn encode_delta_varint<W: Write>(
    posts: &[Posting],
    heads: &[u32],
    mut w: W,
) -> std::io::Result<usize> {
    let mut bytes = 0usize;
    let mut runs = heads
        .iter()
        .copied()
        .map(|x| x as usize)
        .collect::<Vec<_>>();
    runs.push(posts.len());
    for r in 0..(runs.len() - 1) {
        let start = runs[r];
        let end = runs[r + 1];
        if start >= end {
            continue;
        }
        // sort per-run for better delta compression
        let mut buf = posts[start..end].to_vec();
        buf.sort_by_key(|p| ((p.read_id as u64) << 32) | (p.pos as u64));

        let mut prev_id: u32 = 0;
        let mut prev_pos: u32 = 0;
        let mut first = true;
        for p in &buf {
            if first {
                // store absolute first
                bytes += write_varu32(p.read_id, &mut w)?;
                bytes += write_varu32(p.pos, &mut w)?;
                prev_id = p.read_id;
                prev_pos = p.pos;
                first = false;
            } else {
                let d_id = p.read_id.wrapping_sub(prev_id);
                let d_pos = p.pos.wrapping_sub(prev_pos);
                bytes += write_varu32(d_id, &mut w)?;
                bytes += write_varu32(d_pos, &mut w)?;
                prev_id = p.read_id;
                prev_pos = p.pos;
            }
        }
    }
    Ok(bytes)
}

/// Decode postings for an entire bucket from compressed bytes.
///
/// Note: we decode sequentially until EOF of the bucket region; `heads` is
/// present for symmetry but not explicitly used here to split runs.
pub fn decode_delta_varint(
    bytes: &[u8],
    _heads: &[u32],
    out: &mut Vec<Posting>,
) -> std::io::Result<()> {
    use std::io::Cursor;
    let mut rdr = Cursor::new(bytes);

    let mut prev_id = 0u32;
    let mut prev_pos = 0u32;
    let mut first_in_run = true;

    loop {
        let id = match read_varu32(&mut rdr) {
            Ok(v) => v,
            Err(_) => break,
        };
        let pos = match read_varu32(&mut rdr) {
            Ok(v) => v,
            Err(_) => break,
        };

        if first_in_run {
            // absolute
            out.push(Posting { read_id: id, pos });
            prev_id = id;
            prev_pos = pos;
            first_in_run = false;
        } else {
            let rid = prev_id.wrapping_add(id);
            let p = prev_pos.wrapping_add(pos);
            out.push(Posting {
                read_id: rid,
                pos: p,
            });
            prev_id = rid;
            prev_pos = p;
        }
    }
    Ok(())
}

/// `(key,pos)` sample pair for PLR sampling table.
#[repr(C)]
#[derive(Copy, Clone, Default, Pod, Zeroable)]
struct SamplePair {
    key: u64,
    pos: u64,
}
