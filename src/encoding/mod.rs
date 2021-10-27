//! Mod define encoding trait and type implementing this trait

/* project use */

/* mod declaration */
pub mod naive;
pub mod xor10;

/* public use */
pub use naive::Naive;
pub use xor10::Xor10;

/// Trait use by Kmer struct, to convert DNA encode on 8 bits to 2 bits encoding, inverte this operation and perform a reverse complement on the 2 bits encoding
pub trait Encoding<P, const B: usize> {
    /// Convert a DNA sequence, encode with 8 bits per nucleotide in a DNA sequence encode on 2 bits per nucleotide
    fn encode(&self, seq: &[u8]) -> [P; B];

    /// Convert a DNA sequence, encode on 2 bits per nucleotide in a DNA sequence on 8 bits per nucleotide
    fn decode(&self, array: [P; B]) -> Vec<u8>;

    /// Perform a reverse complement on a DNA sequence encode on 2 bits per nucleotide
    fn rev_comp<const K: usize>(&self, array: [P; B]) -> [P; B];
}
