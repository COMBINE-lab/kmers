//! Many utils function

/// Compute the number of possible kmer for a value of k
pub const fn kmer_space(k: u32) -> usize {
    2_usize.pow(k * 2)
}

/// Compute the number of possible canonical kmer for a value of k
pub const fn canonical_space(k: u32) -> usize {
    if k % 2 == 1 {
        2_usize.pow(k * 2) / 2
    } else {
        2_usize.pow(k * 2) / 2 - (k * 2) as usize
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::kmer_space;

    use super::*;

    #[test]
    fn kmer_space_size() {
        assert_eq!(kmer_space(1), 4);
        assert_eq!(kmer_space(2), 16);
        assert_eq!(kmer_space(3), 64);
        assert_eq!(kmer_space(4), 256);
        assert_eq!(kmer_space(5), 1024);
        assert_eq!(kmer_space(6), 4096);
        assert_eq!(kmer_space(7), 16384);
        assert_eq!(kmer_space(8), 65536);
        assert_eq!(kmer_space(9), 262144);
        assert_eq!(kmer_space(10), 1048576);
    }

    #[test]
    fn canonical_space_size() {
        assert_eq!(canonical_space(0), 0);
        assert_eq!(canonical_space(1), 2);
        assert_eq!(canonical_space(2), 4);
        assert_eq!(canonical_space(3), 32);
        assert_eq!(canonical_space(4), 120);
        assert_eq!(canonical_space(5), 512);
        assert_eq!(canonical_space(6), 2036);
        assert_eq!(canonical_space(7), 8192);
        assert_eq!(canonical_space(8), 32752);
        assert_eq!(canonical_space(9), 131072);
        assert_eq!(canonical_space(10), 524268);
    }
}
