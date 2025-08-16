use kira_cdh_compat_kmer_indexer::encode::*;

#[test]
fn test_encode_revcomp_canonical() {
    let s = b"AC";
    let k = 2;
    let code = encode_kmer(s).unwrap(); // LSB
    assert_eq!(code, 0b0001);

    let rc = revcomp(code, k);
    assert_eq!(rc, 0b1011);

    let can = canonical(code, k); // MSB-aligned
    assert_eq!(can, (0b0001u64) << (64 - 2 * k as u64));
}
