use super::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Eq, PartialEq, Default, Debug, Clone, Ord, PartialOrd)]
pub struct Kmer {
    pub k: u8,
    pub(crate) data: u64,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone, Copy)]
pub enum Orientation {
    IsCanonical,
    NotCanononical,
}

const BASE_TABLE: [char; 4] = ['a', 'c', 'g', 't'];
// const RC_TABLE: [u64; 4] = [0x11, 0x10, 0x01, 0x00];

// A bitmask that masks out the topmost 64-pos bits
// Taking the bitwise `&` with this bitmask will return
// the number setting all bits above the `pos` bit to 0.
const fn bitmask(pos: u64) -> u64 {
    (1 << pos) - 1
}
impl Kmer {
    pub fn len(&self) -> usize {
        self.k as usize
    }

    pub fn is_empty(&self) -> bool {
        // Satisfying clippy/rust best practices
        // not sure why we need this right now though...
        self.k == 0
    }

    pub fn from_u64(data: u64, k: u8) -> Self {
        let data = data & MASK_TABLE[k as usize];
        Kmer { data, k }
    }

    pub fn into_u64(&self) -> u64 {
        // easier than using .clone() and From
        self.data
    }

    pub fn is_canonical(&self) -> bool {
        let rc = self.to_reverse_complement();
        *self <= rc
    }

    pub fn orientation(&self) -> Orientation {
        if self.is_canonical() {
            Orientation::IsCanonical
        } else {
            Orientation::NotCanononical
        }
    }

    pub fn to_canonical(&self) -> Self {
        if self.is_canonical() {
            self.clone()
        } else {
            self.to_reverse_complement()
        }
    }

    #[inline]
    pub fn prepend_base_u8(&mut self, c: u8) -> Base {
        let r = (self.data >> (2 * self.k - 2)) & 0x03;
        self.data = MASK_TABLE[self.k as usize] & ((self.data << 2) | encode_binary_u8(c));
        r
    }

    #[inline]
    pub fn append_base_u8(&mut self, c: u8) -> Base {
        let r = self.data & 0x03;
        self.data = (self.data >> 2) | (encode_binary_u8(c) << (2 * self.k - 2));
        r
    }

    #[inline]
    pub fn prepend_base(&mut self, c: Base) -> Base {
        let r = (self.data >> (2 * self.k - 2)) & 0x03;
        self.data = MASK_TABLE[self.k as usize] & ((self.data << 2) | c);
        r
    }

    #[inline]
    pub fn append_base(&mut self, c: Base) -> Base {
        let r = self.data & 0x03;
        self.data = (self.data >> 2) | (c << (2 * self.k - 2));
        r
    }

    /*
         * for now, use ugly but fast version below
        pub fn to_reverse_complement(&self) -> Self {
        let mut old = self.data;
        let mut new = 0u64;
        for _ in 0..self.k {
        let b = old & 3;
        let cb = complement_base(b);
        new <<= 2;
        new |= cb;
        old >>= 2;
    }
        Self {
        data: new,
        k: self.k,
    }
    }
         */

    // adapted from https://www.biostars.org/p/113640/
    pub fn to_reverse_complement(&self) -> Self {
        let mut res = !(self.data);
        res = (res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2;
        res = (res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4;
        res = (res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8;
        res = (res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16;
        res = (res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32;

        Self {
            data: res >> (2 * (32 - self.k)),
            k: self.k,
        }
    }

    pub fn get_reverse_complement_word(w: u64, k: u8) -> u64 {
        let mut res = !w;
        res = (res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2;
        res = (res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4;
        res = (res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8;
        res = (res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16;
        res = (res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32;

        res >> (2 * (32 - k))
    }
}

impl From<Kmer> for String {
    fn from(kmer: Kmer) -> Self {
        let mut s = String::with_capacity(kmer.k as usize);
        let mut w = kmer.data;
        for _ in 0..kmer.k {
            let c_i = w & 3u64; //mask lowerorder bits
            s.push(BASE_TABLE[c_i as usize]);
            w >>= 2;
        }
        s
    }
}

impl From<String> for Kmer {
    fn from(s: String) -> Self {
        if s.len() > 32 {
            panic!("kmers longer than 32 bases not supported");
        }

        let k = s.len() as u8;

        let mut w = 0_u64;
        // read sequence "left to right" from "lower to higher" order bits
        for c in s.chars().rev() {
            w <<= 2;
            w |= encode_binary(c);
        }
        let data = w;
        Kmer { data, k }
    }
}

impl From<&[u8]> for Kmer {
    fn from(s: &[u8]) -> Self {
        if s.len() > 32 {
            panic!("kmers longer than 32 bases not supported");
        }

        let k = s.len() as u8;

        let mut w = 0_u64;
        // read sequence "left to right" from "lower to higher" order bits
        for c in s.iter().rev() {
            w <<= 2;
            w |= encode_binary(*c as char);
        }
        let data = w;
        Kmer { data, k }
    }
}

impl From<&str> for Kmer {
    fn from(s: &str) -> Self {
        Self::from(String::from(s))
    }
}

impl From<Kmer> for u64 {
    fn from(kmer: Kmer) -> Self {
        kmer.data
    }
}

impl std::fmt::Display for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let msg: String = self.clone().into();
        write!(f, "{}", msg)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    //#[quickcheck]
    fn rc_identity(word: u64) -> bool {
        let km = Kmer::from_u64(word, 31);
        km == km.to_reverse_complement().to_reverse_complement()
    }

    //#[quickcheck]
    fn to_canonical_is_canonical(word: u64) -> bool {
        let km = Kmer::from_u64(word, 31);
        km.to_canonical().is_canonical()
    }

    #[test]
    fn test_into_canon() {
        let seq1 = Kmer::from("taa");
        let seq2 = Kmer::from("tta");
        assert_eq!(seq1.to_canonical(), seq1);
        assert_eq!(seq2.to_canonical(), seq1);

        let seq1 = Kmer::from("atc");
        let seq2 = Kmer::from("gat");
        assert_eq!(seq1.to_canonical(), seq1);
        assert_eq!(seq2.to_canonical(), seq1);

        let not_canon = Kmer::from("gatacataggatgg");
        let rc = Kmer::from("gatacataggatgg").to_reverse_complement();

        assert_eq!(rc, not_canon.to_canonical());

        let canon = Kmer::from("agatacataggatgg");
        assert_eq!(canon, canon.to_canonical());
    }

    #[test]
    fn test_is_canon() {
        assert!(Kmer::from("agatacataggatgg").is_canonical());
        assert!(!Kmer::from("gatacataggatgg").is_canonical());
    }

    #[test]
    fn test_ord() {
        assert!(Kmer::from("tcc") < Kmer::from("cct"));
    }

    #[test]
    fn test_append() {
        let mut k1 = Kmer::from("att");
        let k2 = Kmer::from("ttc");

        let shift_off = k1.append_base_u8(b'c');
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, A);

        let mut k1 = Kmer::from("ttcga");
        let k2 = Kmer::from("tcgag");

        let shift_off = k1.append_base_u8(b'g');
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, T);

        let mut k1 = Kmer::from("att");
        let k2 = Kmer::from("ttc");

        let shift_off = k1.append_base(encode_binary_u8(b'c'));
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, A);

        let mut k1 = Kmer::from("ttcga");
        let k2 = Kmer::from("tcgag");

        let shift_off = k1.append_base(encode_binary_u8(b'g'));
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, T);
    }

    #[test]
    fn test_prepend() {
        let mut k1 = Kmer::from("att");
        let k2 = Kmer::from("cat");

        let shift_off = k1.prepend_base_u8(b'c');
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, T);

        let mut k1 = Kmer::from("ttcga");
        let k2 = Kmer::from("gttcg");

        let shift_off = k1.prepend_base_u8(b'g');
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, A);

        let mut k1 = Kmer::from("att");
        let k2 = Kmer::from("cat");

        let shift_off = k1.prepend_base(encode_binary_u8(b'c'));
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, T);

        let mut k1 = Kmer::from("ttcga");
        let k2 = Kmer::from("gttcg");

        let shift_off = k1.prepend_base(encode_binary_u8(b'g'));
        assert_eq!(k1, k2, "{} {}", k1, k2);
        assert_eq!(shift_off, A);
    }

    #[test]
    fn test_rc() {
        let rc = Kmer { k: 1, data: 0 };
        let rc = rc.to_reverse_complement();
        let checked = Kmer::from("t");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("a").to_reverse_complement();
        let checked = Kmer::from("t");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("aaa").to_reverse_complement();
        let checked = Kmer::from("ttt");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("ttt").to_reverse_complement();
        let checked = Kmer::from("aaa");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("ta").to_reverse_complement();
        let checked = Kmer::from("ta");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("ccg").to_reverse_complement();
        let checked = Kmer::from("cgg");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("aat").to_reverse_complement();
        let checked = Kmer::from("att");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("aat").to_reverse_complement();
        let checked = Kmer::from("att");
        assert_eq!(rc, checked, "{} {}", rc, checked);

        let rc = Kmer::from("gatacataggatgg").to_reverse_complement();
        let checked = Kmer::from("ccatcctatgtatc");
        assert_eq!(rc, checked, "{} {}", rc, checked);
    }

    #[test]
    fn str_repr() {
        let seq = String::from("catagatacat");
        let seq_: String = Kmer::from("catagatacat").into();
        assert_eq!(seq, seq_);
    }

    #[test]
    fn bin_repr() {
        let aaa = 0b000000;
        let aac = 0b010000;
        let acc = 0b010100;
        let ccc = 0b010101;
        let aaa_: u64 = Kmer::from("aaa").into();
        let aac_: u64 = Kmer::from("aac").into();
        let acc_: u64 = Kmer::from("acc").into();
        let ccc_: u64 = Kmer::from("ccc").into();

        assert_eq!(aaa, aaa_);
        assert_eq!(aac, aac_);
        assert_eq!(acc, acc_);
        assert_eq!(ccc, ccc_);
    }

    #[test]
    fn aaa() {
        let x = Kmer::from("aaa");
        assert_eq!(x, Kmer::from_u64(0, 3));

        assert_eq!(x.data, 0);
        assert_eq!(x.k, 3);

        let repr: u64 = x.into();
        assert_eq!(repr, 0);

        for k in 1..33 {
            let x = Kmer::from("A".repeat(k));
            assert_eq!(x.data, 0);
            assert_eq!(x.k, k as u8);
        }
    }

    #[test]
    fn test_eq() {
        assert_eq!(Kmer::from("aaa"), Kmer::from("AAA"));
        assert_eq!(Kmer::from("aCa"), Kmer::from("AcA"));

        assert_ne!(Kmer::from("a"), Kmer::from("aa"));
    }

    #[test]
    #[should_panic]
    fn too_long() {
        let _ = Kmer::from("a".repeat(33));
    }

    #[test]
    fn not_too_long() {
        let _ = Kmer::from("a".repeat(32));
    }

    #[test]
    fn test_encode_binary() {
        assert_eq!(encode_binary('A'), A);
        assert_eq!(encode_binary('a'), A);
        assert_eq!(encode_binary('C'), C);
        assert_eq!(encode_binary('c'), C);
        assert_eq!(encode_binary('G'), G);
        assert_eq!(encode_binary('g'), G);
        assert_eq!(encode_binary('T'), T);
        assert_eq!(encode_binary('t'), T);
    }

    #[test]
    #[should_panic]
    fn encode_panics() {
        encode_binary('N');
    }

    #[test]
    fn test_complement_base() {
        assert_eq!(complement_base(A), T);
        assert_eq!(complement_base(T), A);
        assert_eq!(complement_base(C), G);
        assert_eq!(complement_base(G), C);
    }

    #[test]
    fn test_is_valid_nuc() {
        assert!(is_valid_nuc(0));
        assert!(is_valid_nuc(0));
        assert!(is_valid_nuc(0));
        assert!(is_valid_nuc(0));

        assert!(is_valid_nuc(A));
        assert!(is_valid_nuc(C));
        assert!(is_valid_nuc(G));
        assert!(is_valid_nuc(T));

        assert!(!is_valid_nuc(5));
        assert!(!is_valid_nuc(3112));
    }
}

// table that contains bit patterns to mask out the top bits of a word.
// The table is such that maskTable[k] will mask out the top (64 - 2*k) bits of
// the word.
const MASK_TABLE: [u64; 33] = [
    bitmask(0),
    bitmask(2),
    bitmask(4),
    bitmask(6),
    bitmask(8),
    bitmask(10),
    bitmask(12),
    bitmask(14),
    bitmask(16),
    bitmask(18),
    bitmask(20),
    bitmask(22),
    bitmask(24),
    bitmask(26),
    bitmask(28),
    bitmask(30),
    bitmask(32),
    bitmask(34),
    bitmask(36),
    bitmask(38),
    bitmask(40),
    bitmask(42),
    bitmask(44),
    bitmask(46),
    bitmask(48),
    bitmask(50),
    bitmask(52),
    bitmask(54),
    bitmask(56),
    bitmask(58),
    bitmask(60),
    bitmask(62),
    0,
];
