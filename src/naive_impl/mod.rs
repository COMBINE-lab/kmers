mod canonical_kmer;
pub mod canonical_kmer_iterator;
mod kmer;
mod seq_vector;

// re-exports
pub use canonical_kmer::{CanonicalKmer, MatchType};
pub use canonical_kmer_iterator::CanonicalKmerIterator;
pub use kmer::Kmer;

pub use prelude::Base;
pub use prelude::{A, C, G, T};

pub mod prelude {
    pub type Base = u64;
    pub const A: Base = 0;
    pub const C: Base = 1;
    pub const G: Base = 2;
    pub const T: Base = 3;

    #[inline]
    pub fn encode_binary(c: char) -> Base {
        // might have to play some tricks for lookup in a const
        // array at some point
        match c {
            'A' | 'a' => A,
            'C' | 'c' => C,
            'G' | 'g' => G,
            'T' | 't' => T,
            _ => panic!("cannot decode {c} into 2 bit encoding"),
        }
    }

    #[inline]
    pub fn encode_binary_u8(c: u8) -> Base {
        // might have to play some tricks for lookup in a const
        // array at some point
        match c {
            b'A' | b'a' => A,
            b'C' | b'c' => C,
            b'G' | b'g' => G,
            b'T' | b't' => T,
            _ => u64::MAX,
        }
    }

    #[allow(dead_code)]
    #[inline]
    pub fn encode_complement_binary_u8(c: u8) -> Base {
        // might have to play some tricks for lookup in a const
        // array at some point
        match c {
            b'A' | b'a' => T,
            b'C' | b'c' => G,
            b'G' | b'g' => C,
            b'T' | b't' => A,
            _ => u64::MAX,
        }
    }

    #[allow(dead_code)]
    #[inline]
    pub fn encode_complement_binary(c: char) -> Base {
        // might have to play some tricks for lookup in a const
        // array at some point
        match c {
            'A' | 'a' => T,
            'C' | 'c' => G,
            'G' | 'g' => C,
            'T' | 't' => A,
            _ => panic!("cannot decode {c} into 2 bit encoding"),
        }
    }

    #[inline]
    pub fn complement_base(b: Base) -> Base {
        // this is cool
        3 - b
    }

    #[allow(dead_code)]
    pub fn is_valid_nuc(b: Base) -> bool {
        b < 4
    }
}
