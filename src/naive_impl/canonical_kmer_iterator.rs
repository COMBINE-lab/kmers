// This code is derived and modified from
// https://github.com/COMBINE-lab/pufferfish/blob/develop/include/CanonicalKmerIterator.hpp
// which is itself derived and modified from Pall Melsted's
// KmerIterator class (part of bfgraph) :
// https://github.com/pmelsted/bfgraph/blob/master/src/KmerIterator.cpp

use super::prelude::*;
use super::CanonicalKmer;
use serde::{Deserialize, Serialize};

// holds what is essentially a pair of
// km: the canonical k-mer on the read
// pos: the offset on the read where this k-mer starts
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CanonicalKmerPos {
    pub km: CanonicalKmer,
    pub pos: i32,
}

impl CanonicalKmerPos {
    fn new(k: u8) -> Self {
        Self {
            km: CanonicalKmer::blank_of_size(k),
            pos: -1i32,
        }
    }
}

// A CanonicalKmerIterator holds a
// shared reference to an underlying [u8] slice `seq`.
// It is capable of iterating over this sequence (skipping invalid
// k-mers, e.g. k-mers containing `N`), and producing
// a `CanonicalKmerPos` struct for all valid k-mers in `seq`.
#[derive(Debug, Clone)]
pub struct CanonicalKmerIterator<'a> {
    seq: &'a [u8],
    value_pair: CanonicalKmerPos,
    invalid: bool,
    last_invalid: i32,
    k: i32,
}

impl<'slice> CanonicalKmerIterator<'slice> {
    #[inline]
    fn find_next(&mut self, ii: i32, jj: i32) {
        let mut i = ii + 1;
        let j = jj + 1;

        let seq_len = self.seq.len() as i32;

        // l is the last nucleotide in the k-mer we are
        // currently building
        for l in j..seq_len {
            // get the code for the last nucleotide, save it as b
            let b = encode_binary_u8(self.seq[l as usize]);

            // c is an invalid code if >= 4
            if b < 4 {
                self.value_pair.km.append_base(b);
                if (l - self.last_invalid) >= self.k {
                    self.value_pair.pos = i;
                    return;
                }
            } else {
                // this k-mer is clearly not valid
                // if b is not a valid code, then l is the last invalid position
                self.last_invalid = l;
                i = l + 1;
            }
        }

        self.invalid = true;
    }

    pub fn from_u8_slice(s: &'slice [u8], k: u8) -> CanonicalKmerIterator {
        let mut r = Self {
            seq: s,
            value_pair: CanonicalKmerPos::new(k),
            invalid: false,
            last_invalid: -1i32,
            k: k as i32,
        };

        r.find_next(-1, -1);
        r
    }

    // returns true if this iterator is exhausted
    // (i.e. if there are no more valid k-mers beyond)
    // the current position, and false otherwise.
    #[inline]
    pub fn exhausted(&self) -> bool {
        self.invalid
    }

    #[inline]
    pub fn inc(&mut self) -> bool {
        let lpos = self.value_pair.pos + self.k;
        self.invalid = self.invalid || (lpos >= self.seq.len() as i32);
        if !self.invalid {
            self.find_next(self.value_pair.pos, lpos - 1);
        }
        !self.invalid
    }

    #[inline]
    pub fn inc_by(&mut self, mut count: usize) -> bool {
        let mut v = !self.invalid;
        while count > 0 && v {
            v = self.inc();
            count -= 1;
        }
        v
    }

    #[inline]
    pub fn get(&self) -> &CanonicalKmerPos {
        &self.value_pair
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter_init() {
        let r = b"TTTTGGCCATTTTTCCTGTTCTTCAAGAAAACAGGAGATAACTAGAAGGACTAGAGAATGGGGCTGCCAGAACTAGTGGGAAGCTCCCTAGAAATGGTGACATCGCCCACCAAACAGACC";

        let k = 31u8;
        let fk = CanonicalKmer::from(&r[0..31]);

        let ck_iter = CanonicalKmerIterator::from_u8_slice(&r[..], k);

        assert_eq!(fk, ck_iter.get().km);
        assert_eq!(0, ck_iter.get().pos);
    }

    #[test]
    fn test_iter_inc() {
        let r = b"TTTTGGCCATTTTTCCTGTTCTTCAAGAAAACAGGAGATAACTAGAAGGACTAGAGAATGGGGCTGCCAGAACTAGTGGGAAGCTCCCTAGAAATGGTGACATCGCCCACCAAACAGACC";

        let k = 31u8;
        let fk = CanonicalKmer::from(&r[1..32]);

        let mut ck_iter = CanonicalKmerIterator::from_u8_slice(&r[..], k);
        ck_iter.inc();

        assert_eq!(fk, ck_iter.get().km);
        assert_eq!(1, ck_iter.get().pos);
    }

    #[test]
    fn test_iter_inc_by() {
        let r = b"TTTTGGCCATTTTTCCTGTTCTTCAAGAAAACAGGAGATAACTAGAAGGACTAGAGAATGGGGCTGCCAGAACTAGTGGGAAGCTCCCTAGAAATGGTGACATCGCCCACCAAACAGACC";

        let k = 31u8;
        let fk = CanonicalKmer::from(&r[10..41]);

        let mut ck_iter = CanonicalKmerIterator::from_u8_slice(&r[..], k);
        ck_iter.inc_by(10);

        assert_eq!(fk, ck_iter.get().km);
        assert_eq!(10, ck_iter.get().pos);
    }

    #[test]
    fn test_iter_init_invalid() {
        let r = b"TTTTNGGCCATTTTTCCTGTTCTTCAAGAAAACAGGAGATAACTAGAAGGACTAGAGAATGGGGCTGCCAGAACTAGTGGGAAGCTCCCTAGAAATGGTGACATCGCCCACCAAACAGACC";

        let k = 31u8;
        let fk = CanonicalKmer::from(&r[5..36]);

        let ck_iter = CanonicalKmerIterator::from_u8_slice(&r[..], k);

        assert_eq!(fk, ck_iter.get().km);
        assert_eq!(5, ck_iter.get().pos);
    }

    #[test]
    fn test_iter_inc_by_invalid() {
        let r = b"TTTTGGCCATTTTTCCTGTTCTTCAAGAAAACAGGNAGATAACTAGAAGGACTAGAGAATGGGGCTGCCAGAACTAGTGGGAAGCTCCCTAGAAATGGTGACATCGCCCACCAAACAGACC";

        let k = 31u8;
        let fk = CanonicalKmer::from(&r[36..67]);

        let mut ck_iter = CanonicalKmerIterator::from_u8_slice(&r[..], k);
        ck_iter.inc_by(5);

        assert_eq!(fk, ck_iter.get().km);
        assert_eq!(36, ck_iter.get().pos);
    }

    #[test]
    fn test_exhausted_works() {
        let r = b"TTTTGGCCATTTTTCCTGTTCTTCAAGAAAACAGGAGATAACTAGAAGGACTAGAGAATGGGGCTGCCAGAACTAGTGGGAAGCTCCCTAGAAATGGTGACATCGCCCACCAAACAGACC";
        let sl = r.len();
        let k = 31u8;
        let mut ck_iter = CanonicalKmerIterator::from_u8_slice(&r[..], k);
        ck_iter.inc_by(20);

        assert!(!ck_iter.exhausted());

        ck_iter.inc_by(sl - 20);
        assert!(ck_iter.exhausted());

        ck_iter.inc();
        assert!(ck_iter.exhausted());
    }
}
