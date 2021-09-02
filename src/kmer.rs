//! Mod define generic struct to store kmer

/* crate use */
use bit_field::BitArray;

/// Struct to store and use kmer
pub struct Kmer<P, const K: usize, const B: usize> {
    array: [P; B],
}

impl<P, const K: usize, const B: usize> Kmer<P, K, B>
where
    P: Copy + bit_field::BitField,
{
    /// construct a new empty k-mer fill with zero by default
    pub fn new() -> Self {
        Self {
            array: [unsafe { std::mem::zeroed() }; B],
        }
    }

    /// returns the value of k for this k-mer
    pub fn k(&self) -> usize {
        K as usize
    }

    /// returns the number of bytes used for the storage of this k-mer
    pub fn num_bytes(&self) -> usize {
        std::mem::size_of::<P>() * word_for_k::<P, K>()
    }

    /// get the niest nucleotide
    pub fn get(&self, index: usize) -> P {
        self.array.get_bits(index..index + 1)
    }
}

impl<P, const K: usize, const B: usize> std::default::Default for Kmer<P, K, B>
where
    P: Copy + bit_field::BitField,
{
    fn default() -> Self {
        Self::new()
    }
}

/// compute the number of words required to store a kmer of length k
pub const fn word_for_k<P, const K: usize>() -> usize {
    (std::mem::size_of::<P>() * 8 / 2 + K - 1) / (std::mem::size_of::<P>() * 8 / 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn choose_number_of_word() {
        assert_eq!(word_for_k::<u8, 1>(), 1);
        assert_eq!(word_for_k::<u8, 4>(), 1);
        assert_eq!(word_for_k::<u8, 5>(), 2);

        assert_eq!(word_for_k::<u16, 1>(), 1);
        assert_eq!(word_for_k::<u16, 8>(), 1);
        assert_eq!(word_for_k::<u16, 9>(), 2);

        assert_eq!(word_for_k::<u32, 1>(), 1);
        assert_eq!(word_for_k::<u32, 16>(), 1);
        assert_eq!(word_for_k::<u32, 17>(), 2);

        assert_eq!(word_for_k::<u64, 1>(), 1);
        assert_eq!(word_for_k::<u64, 32>(), 1);
        assert_eq!(word_for_k::<u64, 64>(), 2);

        assert_eq!(word_for_k::<u128, 1>(), 1);
        assert_eq!(word_for_k::<u128, 64>(), 1);
        assert_eq!(word_for_k::<u128, 65>(), 2);
    }

    #[test]
    fn u8_k15() {
        let kmer = Kmer::<u8, 15, { word_for_k::<u8, 15>() }>::new();
        assert_eq!(kmer.num_bytes(), 4);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u16_k15() {
        let kmer = Kmer::<u16, 15, { word_for_k::<u16, 15>() }>::new();
        assert_eq!(kmer.num_bytes(), 4);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u32_k15() {
        let kmer = Kmer::<u32, 15, { word_for_k::<u32, 15>() }>::new();
        assert_eq!(kmer.num_bytes(), 4);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u64_k15() {
        let kmer = Kmer::<u64, 15, { word_for_k::<u64, 15>() }>::new();
        assert_eq!(kmer.num_bytes(), 8);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u128_k15() {
        let kmer = Kmer::<u128, 15, { word_for_k::<u128, 15>() }>::new();
        assert_eq!(kmer.num_bytes(), 16);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn get() {
        let kmer = Kmer::<u8, 15, { word_for_k::<u8, 15>() }>::new();

        assert_eq!(kmer.get(2), 0b00)
    }
}
