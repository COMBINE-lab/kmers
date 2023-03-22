use std::hash::BuildHasher;

use serde::{Deserialize, Serialize};
use simple_sds::int_vector::IntVector;
use simple_sds::ops::Vector;
use simple_sds::raw_vector::{AccessRaw, PushRaw, RawVector};

use crate::naive_impl::Kmer;
use simple_sds::serde_compat;

use self::minimizers::SeqVecMinimizerIter;

pub mod minimizers;

#[allow(non_camel_case_types)]
type km_size_t = usize;

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeqVector {
    #[serde(with = "serde_compat")]
    data: RawVector,
}

#[derive(Clone, Debug)]
pub struct SeqVectorSlice<'a> {
    len: usize,
    start_pos: usize,
    slice: &'a SeqVector,
}

impl SeqVectorSlice<'_> {
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get_kmer(&self, pos: usize, k: km_size_t) -> Kmer {
        let km = self.get_kmer_u64(pos, k);
        Kmer::from_u64(km, k as u8)
    }

    pub fn get_kmer_u64(&self, pos: usize, k: km_size_t) -> u64 {
        assert!(pos < self.len());
        let pos = pos + self.start_pos;
        self.slice.get_kmer_u64(pos, k)
    }

    pub fn get_base(&self, pos: usize) -> u64 {
        self.get_kmer_u64(pos, 1)
    }

    pub fn slice(&self, start: usize, end: usize) -> Self {
        assert!(end <= self.len());
        Self {
            len: end - start,
            start_pos: start,
            slice: self.slice,
        }
    }

    pub fn iter_kmers(&self, k: km_size_t) -> SeqVecKmerIterator {
        SeqVecKmerIterator {
            k,
            len: self.len - k + 1,
            pos: 0,
            seq: self.clone(),
        }
    }

    pub fn iter_minimizers<T: BuildHasher>(
        &self,
        k: km_size_t,
        w: km_size_t,
        build_hasher: T,
    ) -> SeqVecMinimizerIter<T> {
        SeqVecMinimizerIter::new(self.clone(), k, w, build_hasher)
    }
}

impl SeqVector {
    pub fn len(&self) -> usize {
        self.data.len() / 2
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn get_kmer(&self, pos: usize, k: km_size_t) -> Kmer {
        Kmer::from_u64(self.get_kmer_u64(pos, k), k as u8)
    }

    pub fn get_kmer_u64(&self, pos: usize, k: km_size_t) -> u64 {
        assert!(pos < self.len());
        unsafe { self.data.int(pos * 2, k * 2) }
    }

    pub fn get_base(&self, pos: usize) -> u64 {
        self.get_kmer_u64(pos, 1)
    }

    pub fn as_slice(&self) -> SeqVectorSlice<'_> {
        SeqVectorSlice {
            start_pos: 0,
            len: self.len(),
            slice: self,
        }
    }

    pub fn new() -> Self {
        Self {
            data: RawVector::new(),
        }
    }

    pub fn with_len(len: usize) -> Self {
        // initializes with all 0b00s, i.e. As.
        Self {
            data: RawVector::with_len(len * 2, false),
        }
    }

    pub fn slice(&self, start: usize, end: usize) -> SeqVectorSlice {
        self.as_slice().slice(start, end)
    }

    pub fn iter_kmers(&self, k: km_size_t) -> SeqVecKmerIterator {
        SeqVecKmerIterator {
            k,
            len: self.len() - k + 1,
            pos: 0,
            seq: self.as_slice(),
        }
    }

    pub fn iter_minimizers<T: BuildHasher>(
        &self,
        k: km_size_t,
        w: km_size_t,
        build_hasher: T,
    ) -> SeqVecMinimizerIter<T> {
        SeqVecMinimizerIter::new(self.as_slice(), k, w, build_hasher)
    }

    pub fn with_capacity(len: usize) -> Self {
        Self {
            data: RawVector::with_capacity(len * 2),
        }
    }

    pub fn set_chars(&mut self, offset: usize, bytes: &[u8]) {
        assert!(offset + bytes.len() <= self.len());

        let first_word_len = 32 - (offset % 32);
        let first_word_len = usize::min(first_word_len, bytes.len());

        let (first, rest) = bytes.split_at(first_word_len);

        let last_word_len = rest.len() % 32;
        let (rest, last) = rest.split_at(rest.len() - last_word_len);

        let mut offset = offset;

        if !first.is_empty() {
            let first = Kmer::from(first).into_u64();
            // push the first
            unsafe {
                self.data.set_int(offset * 2, first, first_word_len * 2);
            }
            offset += first_word_len;
        }

        // push the rest that is u64 aligned.
        let chunks = rest.chunks(32);
        for chunk in chunks {
            let word = Kmer::from(chunk).into_u64();
            unsafe {
                self.data.set_int(offset * 2, word, chunk.len() * 2);
            }
            offset += chunk.len();
        }

        if !last.is_empty() {
            let last = Kmer::from(last).into_u64();
            // push the first
            unsafe {
                self.data.set_int(offset * 2, last, last_word_len * 2);
            }
        }
    }

    pub fn push_chars(&mut self, bytes: &[u8]) {
        // push chars so that they are u64 aligned

        let first_word_len = 32 - (self.len() % 32);
        let first_word_len = usize::min(first_word_len, bytes.len());

        let (first, rest) = bytes.split_at(first_word_len);

        let last_word_len = rest.len() % 32;
        let (rest, last) = rest.split_at(rest.len() - last_word_len);

        if !first.is_empty() {
            let first = Kmer::from(first).into_u64();
            // push the first
            unsafe {
                self.data.push_int(first, first_word_len * 2);
            }
        }

        // push the rest that is u64 aligned.
        let chunks = rest.chunks(32);
        for chunk in chunks {
            let word = Kmer::from(chunk).into_u64();
            unsafe {
                self.data.push_int(word, chunk.len() * 2);
            }
        }

        if !last.is_empty() {
            let last = Kmer::from(last).into_u64();
            // push the first
            unsafe {
                self.data.push_int(last, last_word_len * 2);
            }
        }
    }
}

impl std::fmt::Display for SeqVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // write!(f, "SeqVector[ {} ]", String::from(self))
        write!(f, "{}", String::from(self))
    }
}

impl From<&SeqVector> for String {
    fn from(data: &SeqVector) -> Self {
        let mut str = String::new();
        let bases = vec!['A', 'C', 'G', 'T'];
        for i in 0..data.len() {
            let base = data.get_base(i);
            let base = bases[base as usize];
            str.push(base);
        }
        str
    }
}

impl std::fmt::Display for SeqVectorSlice<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // write!(f, "SeqVector[ {} ]", String::from(self))
        write!(f, "{}", String::from(self))
    }
}

impl From<&SeqVectorSlice<'_>> for String {
    fn from(data: &SeqVectorSlice<'_>) -> Self {
        let mut str = String::new();
        let bases = vec!['A', 'C', 'G', 'T'];
        for i in 0..data.len() {
            let base = data.get_base(i);
            let base = bases[base as usize];
            str.push(base);
        }
        str
    }
}

impl From<SeqVector> for String {
    fn from(data: SeqVector) -> Self {
        Self::from(&data)
    }
}

impl From<&String> for SeqVector {
    fn from(data: &String) -> Self {
        assert!(data.is_ascii());
        let bytes = data.as_bytes();
        Self::from(bytes)
    }
}

impl From<String> for SeqVector {
    fn from(data: String) -> Self {
        Self::from(&data)
    }
}

impl<const N: usize> From<&[u8; N]> for SeqVector {
    fn from(data: &[u8; N]) -> Self {
        Self::from(data.as_slice())
    }
}

impl From<&[u8]> for SeqVector {
    fn from(data: &[u8]) -> Self {
        let len = data.len() * 2;
        let chunks = data.chunks(32);
        let mut words = Vec::with_capacity(chunks.len());
        for chunk in chunks {
            let word = Kmer::from(chunk);
            words.push(word.into_u64());
        }
        let rv = RawVector::from_parts(len, words);
        Self { data: rv }
    }
}

impl From<RawVector> for SeqVector {
    fn from(data: RawVector) -> Self {
        assert_eq!(data.len() % 2, 0);
        Self { data }
    }
}

impl From<IntVector> for SeqVector {
    fn from(data: IntVector) -> Self {
        assert_eq!(data.width(), 2);
        Self {
            data: RawVector::from(data),
        }
    }
}

pub struct SeqVecKmerIterator<'a> {
    k: km_size_t,
    len: usize,
    pos: usize,
    seq: SeqVectorSlice<'a>,
}

impl<'a> SeqVecKmerIterator<'a> {
    pub fn new(slice: SeqVectorSlice<'a>, k: km_size_t) -> Self {
        Self {
            k,
            len: slice.len() - k + 1,
            pos: 0,
            seq: slice,
        }
    }
}

impl SeqVecKmerIterator<'_> {
    pub fn len(&self) -> usize {
        // warn!{"len() to be changed to ExactSizeIterator::len(...)"}
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Iterator for SeqVecKmerIterator<'_> {
    type Item = Kmer;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos < self.len() {
            let km = self.seq.get_kmer(self.pos, self.k);
            self.pos += 1;
            Some(km)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {

    use super::super::hash::LexHasherState;
    use super::*;

    #[test]
    fn seq_slice_test() {
        let bytes = vec![1u64, 2, 3];
        let iv = IntVector::from(bytes);
        let rv = RawVector::from(iv);
        let sv = SeqVector::from(rv);

        let slice = sv.as_slice();

        assert_eq!(slice.len(), 32 * 3);
        assert_eq!(slice.get_kmer_u64(0, 32), 1);

        let slice = sv.slice(1, 96);
        assert_eq!(slice.get_kmer_u64(0, 32), sv.get_kmer_u64(1, 32));

        let slice = sv.slice(75, 96);
        assert_eq!(slice.get_kmer_u64(0, 7), sv.get_kmer_u64(75, 7));
    }

    #[test]
    fn push_chars() {
        let mut sv = SeqVector::with_capacity(64);
        let first_a30 = "A".repeat(30);
        let last_c40 = "C".repeat(40);
        sv.push_chars(first_a30.as_bytes());

        assert_eq!(sv.to_string(), first_a30);
        assert_eq!(sv.len(), 30);
        sv.push_chars(last_c40.as_bytes());
        assert_eq!(sv.len(), 70);
        assert_eq!(sv.to_string(), first_a30 + &last_c40);
    }

    #[test]
    fn set_chars() {
        let mut sv = SeqVector::with_capacity(64);
        let first_a30 = "A".repeat(30);
        let last_c40 = "C".repeat(40);
        sv.push_chars(first_a30.as_bytes());
        sv.push_chars(last_c40.as_bytes());

        let set_g40 = "G".repeat(40);
        sv.set_chars(5, set_g40.as_bytes());

        assert_eq!(sv.len(), 70);
        assert_eq!(
            sv.to_string(),
            "A".repeat(5) + &"G".repeat(40) + &"C".repeat(25)
        );

        let mut sv = SeqVector::with_len(70);
        assert_eq!(sv.len(), 70);
        assert_eq!(sv.to_string(), "A".repeat(70),);

        sv.set_chars(15, set_g40.as_bytes());
        assert_eq!(
            sv.to_string(),
            "A".repeat(15) + &"G".repeat(40) + &"A".repeat(15)
        );

        let mut sv = SeqVector::with_len(32);

        sv.set_chars(1, "G".repeat(2).as_bytes());
        sv.set_chars(3, "C".repeat(29).as_bytes());

        assert_eq!(
            sv.to_string(),
            "A".repeat(1) + &"G".repeat(2) + &"C".repeat(29)
        );

        sv.set_chars(0, "C".repeat(32).as_bytes());
        assert_eq!(sv.to_string(), "C".repeat(32),);
    }

    #[test]
    fn iter_kmers() {
        let s = b"ACTTGAT";
        let sv = SeqVector::from(s);
        let mers = vec!["act", "ctt", "ttg", "tga", "gat"];

        let kmers: Vec<String> = sv.iter_kmers(3).map(|km| km.to_string()).collect();
        assert_eq!(kmers, mers);

        let kmers: Vec<String> = sv
            .slice(1, sv.len() - 1)
            .iter_kmers(3)
            .map(|km| km.to_string())
            .collect();
        assert_eq!(kmers, mers[1..mers.len() - 1]);
    }

    #[test]
    fn iter_minimizers() {
        let s = b"ACTTGAT";
        let sv = SeqVector::from(s);
        let k = 5;
        let w = 3;
        let build_hasher = LexHasherState::new(w);

        let _mmers = sv.iter_minimizers(k, w, build_hasher);

        let mers = vec!["act", "ctt", "ttg", "tga", "gat"];

        let kmers: Vec<String> = sv.iter_kmers(3).map(|km| km.to_string()).collect();
        assert_eq!(kmers, mers);

        let kmers: Vec<String> = sv
            .slice(1, sv.len() - 1)
            .iter_kmers(3)
            .map(|km| km.to_string())
            .collect();
        assert_eq!(kmers, mers[1..mers.len() - 1]);
    }
}
