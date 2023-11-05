use super::{Base, Kmer};

/// This module can only return one error, the EncodeError
#[derive(Debug, PartialEq, Eq)]
pub struct EncodeError;
type Result<T> = std::result::Result<T, EncodeError>;

impl Kmer {
    /// Failable conversion from byte slice to Kmer
    pub fn from_bytes_checked(s: &[u8]) -> Result<Self> {
        if s.len() > 32 {
            panic!("kmers longer than 32 bases not supported");
        }

        let k = s.len() as u8;

        let mut w = 0_u64;
        // read sequence "left to right" from "lower to higher" order bits
        for c in s.iter().rev() {
            w <<= 2;
            w |= encode_binary_checked(*c as char)?;
        }
        let data = w;

        Ok(Kmer { data, k })
    }
}

/// Failable encoding
pub fn encode_binary_checked(c: char) -> Result<Base> {
    let code = CODES[c as usize];
    if code >= 0 {
        Ok(code as Base)
    } else {
        Err(EncodeError)
    }
}

// see Kmer.hpp
const R: i32 = -1;
const I: i32 = -2;
const O: i32 = -3;
const A: i32 = 0;
const C: i32 = 1;
const G: i32 = 2;
const T: i32 = 3;

const CODES: [i32; 256] = [
    O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O,
    O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
];

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn checked_nuc_encoding() {
        let allowed = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'];
        let expected_err = Err(EncodeError);
        for c in 0..(char::MAX as u8) {
            let c = c as char;
            if !allowed.contains(&c) {
                assert_eq!(expected_err, encode_binary_checked(c));
            }
        }

        let nucs: Vec<Result<u64>> = allowed.iter().map(|c| encode_binary_checked(*c)).collect();
        let codes: Vec<Result<u64>> = [0, 1, 2, 3, 0, 1, 2, 3].iter().map(|i| Ok(*i)).collect();
        assert_eq!(nucs, codes);
    }

    #[test]
    fn checked_kmer_encoding() {
        let bytes = b"ANa";
        let kw = Kmer::from_bytes_checked(bytes);

        assert_eq!(kw, Err(EncodeError));

        let bytes = b"acgt";
        let km = Kmer::from(bytes);
        let kw = Kmer::from_bytes_checked(bytes).unwrap();

        assert_eq!(km, kw);
        assert_eq!(kw.len(), 4);
    }
}
