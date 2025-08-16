use kira_cdh_compat_kmer_indexer::*;
use proptest::prelude::*;

/// Naive scan using rolling forward/rc for baseline.
fn naive_scan(
    seq: &[u8],
    k: usize,
    canonical: bool,
) -> std::collections::HashMap<u64, Vec<(u32, u32)>> {
    use std::collections::HashMap;
    let mut h = HashMap::<u64, Vec<(u32, u32)>>::new();
    if seq.len() < k {
        return h;
    }
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let lut = encode::MAP_LUT;

    let mut fwd = 0u64;
    let mut rc = 0u64;
    let mut len = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        let v = lut[b as usize];
        if v > 3 {
            fwd = 0;
            rc = 0;
            len = 0;
            continue;
        }
        let v = v as u64;
        fwd = ((fwd << 2) | v) & mask;
        rc = (rc >> 2) | ((v ^ 0b11) << (2 * (k - 1)));
        len += 1;
        if len >= k {
            let code_lsb = if canonical { fwd.min(rc) } else { fwd };
            let code_msb = encode::align_msb(code_lsb, k);
            h.entry(code_msb).or_default().push((0, (i + 1 - k) as u32));
        }
    }
    h
}

proptest! {
    // Lightweight property tests for encode/bucketing primitives.
    #[test]
    fn prop_encode_and_bucket(
        k in 1usize..=12,
        canonical in any::<bool>(),
        seq in prop::collection::vec(prop::sample::select(b"ACGTN".to_vec()), 32..256)
    ) {
        let B = (2*k).min(12) as u8;
        let baseline = naive_scan(&seq, k, canonical);
        for (code, _occ) in baseline {
            let b = encode::bucket_id_from_msb(code, B);
            prop_assert!(b < (1usize << B));
        }
    }
}
