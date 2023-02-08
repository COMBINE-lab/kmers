//! Mod define generic struct to store kmer

/* crate use */
use bit_field::BitArray;
use std::u32;

/* project use */
use crate::encoding;

/// Struct to store and use kmer
#[derive(Debug)]
pub struct Kmer<P, const K: usize, const B: usize> {
    array: [P; B],
}

impl<P, const K: usize, const B: usize> Kmer<P, K, B>
where
    P: Copy + bit_field::BitField,
{
    /// construct a new empty k-mer fill with zero by default
    pub fn new<E>(sequence: &[u8], encoder: &E) -> Self
    where
        E: encoding::Encoding<P, B>,
    {
        Self {
            array: encoder.encode(sequence),
        }
    }

    /// construct a new k-mer with a slice
    pub fn with_data(data: [P; B]) -> Self {
        Self { array: data }
    }

    /// returns the value of k for this k-mer
    pub fn k(&self) -> usize {
        K
    }

    /// returns the number of bytes used for the storage of this k-mer
    pub fn num_bytes(&self) -> usize {
        std::mem::size_of::<P>() * word_for_k::<P, K>()
    }

    /// get the niest nucleotide
    pub fn get(&self, index: usize) -> P {
        self.array.get_bits(index * 2..=index * 2 + 1)
    }

    pub fn get_prefix(&self, len: usize) -> P {
        self.array.get_bits(0..=(len * 2))
    }
}

impl<P, const K: usize, const B: usize> std::default::Default for Kmer<P, K, B>
where
    P: Copy + bit_field::BitField,
{
    fn default() -> Self {
        Self {
            array: [unsafe { std::mem::zeroed() }; B],
        }
    }
}

/// compute the number of words required to store a kmer of length k
pub const fn word_for_k<P, const K: usize>() -> usize {
    (std::mem::size_of::<P>() * 8 / 2 + K - 1) / (std::mem::size_of::<P>() * 8 / 2)
}

pub fn bitmer_to_bytes(mer: u64, len_in: usize) -> Vec<u8> {
    let mut new_kmer = mer;
    let len = len_in as u32;
    let mut new_kmer_str = Vec::new();

    // First char is in the lowest order two bits. So just mask then shift 2.
    let bitmask = 0b11;
    for _ in 0..len {
        // let new_char = (new_kmer & bitmask) >> offset;
        let new_char = new_kmer & bitmask;
        new_kmer >>= 2;
        new_kmer_str.push(match new_char {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!("Mathematical impossibility"),
        });
    }
    new_kmer_str
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
        let kmer = Kmer::<u8, 15, { word_for_k::<u8, 15>() }>::default();
        assert_eq!(kmer.num_bytes(), 4);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u16_k15() {
        let kmer = Kmer::<u16, 15, { word_for_k::<u16, 15>() }>::default();
        assert_eq!(kmer.num_bytes(), 4);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u32_k15() {
        let kmer = Kmer::<u32, 15, { word_for_k::<u32, 15>() }>::default();
        assert_eq!(kmer.num_bytes(), 4);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u64_k15() {
        let kmer = Kmer::<u64, 15, { word_for_k::<u64, 15>() }>::default();
        assert_eq!(kmer.num_bytes(), 8);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn u128_k15() {
        let kmer = Kmer::<u128, 15, { word_for_k::<u128, 15>() }>::default();
        assert_eq!(kmer.num_bytes(), 16);
        assert_eq!(kmer.k(), 15);
    }

    #[test]
    fn kmer_with_data() {
        let data: [u8; 1] = [0b11100100];

        let kmer = Kmer::<u8, 4, { word_for_k::<u8, 4>() }>::with_data(data);

        assert_eq!(kmer.get(0), 0b00);
        assert_eq!(kmer.get(1), 0b01);
        assert_eq!(kmer.get(2), 0b10);
        assert_eq!(kmer.get(3), 0b11);
    }

    #[test]
    fn kmer_naive_encoder() {
        let encoder = encoding::Naive::ACTG;
        let kmer = Kmer::<u8, 4, { word_for_k::<u8, 4>() }>::new(b"ACTG", &encoder);

        assert_eq!(kmer.get(0), 0b00);
        assert_eq!(kmer.get(1), 0b01);
        assert_eq!(kmer.get(2), 0b10);
        assert_eq!(kmer.get(3), 0b11);

        let encoder = encoding::Naive::TAGC;
        let kmer = Kmer::<u8, 4, { word_for_k::<u8, 4>() }>::new(b"ACTG", &encoder);

        assert_eq!(kmer.get(0), 0b01);
        assert_eq!(kmer.get(1), 0b11);
        assert_eq!(kmer.get(2), 0b00);
        assert_eq!(kmer.get(3), 0b10);
    }

    #[test]
    fn kmer_prefix() {
        let encoder = encoding::Naive::ACGT;
        let kmer = Kmer::<u64, 31, { word_for_k::<u64, 31>() }>::new(b"GTAC", &encoder);

        let pref: u64 = kmer.get_prefix(4);
        assert_eq!(pref, 0b01001110);

        let s = bitmer_to_bytes(pref, 4);
        assert_eq!(b"GTAC".to_vec(), s);
    }

    #[test]
    fn kmer_to_bytes() {
        let pref = 0b01001110;
        let s = bitmer_to_bytes(pref, 4);
        assert_eq!(b"GTAC".to_vec(), s);
    }
}
