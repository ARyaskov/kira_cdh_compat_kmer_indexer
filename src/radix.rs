//! In-place LSD radix sort for `u64` keys with a paired values array.
//! 8-bit passes, 8 rounds. Stable via counting + prefix sums.

use crate::index::Posting;

/// Sort `keys` in ascending order and permute `vals` accordingly.
/// Temporary buffers are allocated once and reused across passes.
pub fn radix_sort_pairs_u64(keys: &mut [u64], vals: &mut [Posting]) {
    debug_assert_eq!(keys.len(), vals.len());
    let n = keys.len();
    if n <= 1 {
        return;
    }

    // Scratch buffers.
    let mut tmp_keys = vec![0u64; n];
    let mut tmp_vals = vec![Posting::default(); n];

    // For each byte [0..7], perform a counting sort pass.
    for pass in 0..8 {
        let shift = pass * 8;
        let mut counts = [0usize; 256];

        // Count occurrences
        for &k in keys.iter() {
            let byte = ((k >> shift) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Prefix sums -> positions
        let mut sum = 0usize;
        for c in counts.iter_mut() {
            let tmp = *c;
            *c = sum;
            sum += tmp;
        }

        // Scatter to tmp (stable)
        for i in 0..n {
            let k = keys[i];
            let b = ((k >> shift) & 0xFF) as usize;
            let pos = counts[b];
            tmp_keys[pos] = k;
            tmp_vals[pos] = vals[i];
            counts[b] = pos + 1;
        }

        // Swap buffers
        keys.copy_from_slice(&tmp_keys);
        vals.copy_from_slice(&tmp_vals);
    }
}
