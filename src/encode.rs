//! K-mer encoding: 2-bit mapping, rolling reverse complement, canonicalization.
//!
//! Conventions
//! - `encode_kmer` returns **LSB-aligned** code (lower `2k` bits).
//! - `canonical(code_lsb, k)` returns **MSB-aligned** code (upper `2k` bits).
//! - Storage/bucketing use **MSB-aligned** codes.

/// 256-entry LUT: ASCII → 2-bit (A=0, C=1, G=2, T/U=3), 0xFF for ambiguous.
pub static MAP_LUT: [u8; 256] = {
    const X: u8 = 0xFF;
    let mut t = [X; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
};

/// 2-bit encoding via LUT: A=00, C=01, G=10, T=11. `None` if ambiguous.
#[inline]
pub fn map_base(b: u8) -> Option<u8> {
    let v = MAP_LUT[b as usize];
    if v <= 3 { Some(v) } else { None }
}

/// Encode a k-mer window to **LSB-aligned** `u64`. None if `k>32` or ambiguous.
#[inline]
pub fn encode_kmer(window: &[u8]) -> Option<u64> {
    let k = window.len();
    if k == 0 || k > 32 {
        return None;
    }
    let mut code: u64 = 0;
    for &b in window {
        let v = map_base(b)? as u64;
        code = (code << 2) | v;
    }
    Some(code)
}

/// Reverse-complement an **LSB-aligned** code (lower `2k` bits used).
#[inline]
pub fn revcomp(code_lsb: u64, k: usize) -> u64 {
    debug_assert!(k <= 32);
    let mut rc: u64 = 0;
    for i in 0..k {
        let base = (code_lsb >> (i * 2)) & 0b11;
        let comp = base ^ 0b11;
        let shift = (k - 1 - i) * 2;
        rc |= comp << shift;
    }
    rc
}

/// Canonicalize an **LSB-aligned** code and return **MSB-aligned** code.
#[inline]
pub fn canonical(code_lsb: u64, k: usize) -> u64 {
    let rc = revcomp(code_lsb, k);
    let c = if code_lsb <= rc { code_lsb } else { rc };
    c << (64 - 2 * k as u64)
}

/// Align an LSB-aligned code to MSB-aligned form (upper `2k` bits used).
#[inline]
pub fn align_msb(code_lsb: u64, k: usize) -> u64 {
    code_lsb << (64 - 2 * k as u64)
}

/// Extract the top-B bits (bucket id) from an MSB-aligned code.
#[inline]
pub fn bucket_id_from_msb(code_msb: u64, b: u8) -> usize {
    let b = b.min(63);
    ((code_msb >> (64 - b as u64)) & ((1u64 << b) - 1)) as usize
}
