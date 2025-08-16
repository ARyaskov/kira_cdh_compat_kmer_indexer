//! Lightweight two-level sampling index (PLR-like) to accelerate locate.
//!
//! We store every `stride`-th (key, pos) pair in ascending key order.
//! Query: upper_bound in samples, then binary search in [pos..pos+stride] range.

#[derive(Clone)]
pub struct PlrSamples {
    pub stride: usize,
    pub keys: Vec<u64>, // MSB-aligned keys (monotone increasing)
    pub pos: Vec<u64>,  // positions in the bucket's codes[]
}

impl PlrSamples {
    pub fn build(codes: &[u64], stride: usize) -> Option<Self> {
        if codes.len() < stride {
            return None;
        }
        let mut keys = Vec::with_capacity((codes.len() + stride - 1) / stride);
        let mut pos = Vec::with_capacity(keys.capacity());
        let mut i = 0usize;
        while i < codes.len() {
            keys.push(codes[i]);
            pos.push(i as u64);
            i = i.saturating_add(stride);
        }
        Some(Self { stride, keys, pos })
    }

    #[inline]
    pub fn range_hint(&self, key: u64, max_len: usize) -> (usize, usize) {
        // upper_bound in keys
        let mut lo = 0usize;
        let mut hi = self.keys.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.keys[mid] <= key {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        let idx = lo.saturating_sub(1);
        let start = if self.keys.is_empty() {
            0
        } else {
            self.pos[idx] as usize
        };
        let end = (start + self.stride + 64).min(max_len); // small guard
        (start, end)
    }
}
