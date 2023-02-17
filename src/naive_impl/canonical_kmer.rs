use super::prelude::*;
use super::Kmer;

use serde::{Deserialize, Serialize};
use std::convert::From;

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone, Copy)]
pub enum MatchType {
    NoMatch,
    IdentityMatch,
    TwinMatch,
}

#[derive(Eq, PartialEq, Default, Debug, Clone, Ord, PartialOrd)]
pub struct CanonicalKmer {
    fw: Kmer,
    rc: Kmer,
}

impl CanonicalKmer {
    #[inline]
    pub fn blank_of_size(k: u8) -> Self {
        let data = 0u64;
        let rc_data = u64::MAX;
        Self {
            fw: Kmer { data, k },
            rc: Kmer { data: rc_data, k },
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.fw.is_empty()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.fw.len()
    }

    #[inline]
    pub fn from_u64(data: u64, k: u8) -> Self {
        let fw = Kmer::from_u64(data, k);
        let rc = fw.to_reverse_complement();
        // let rd = Kmer::get_reverse_complement_word(data, k);
        Self {
            fw,
            rc, // fw: Kmer { data, k },
                // rc: Kmer { data: rd, k },
        }
    }

    // #[inline]
    // pub fn from_kmer(km: Kmer) -> Self {
    //     Self {
    //         rc: km.to_reverse_complement(),
    //         fw: km,
    //     }
    // }

    #[inline]
    pub fn swap(&mut self) {
        std::mem::swap(&mut self.fw.data, &mut self.rc.data)
    }

    #[inline]
    pub fn is_fw_canonical(&self) -> bool {
        self.fw.data < self.rc.data
    }

    #[inline]
    pub fn append_base_u8(&mut self, c: u8) -> Base {
        let b = encode_binary_u8(c);
        let cb = complement_base(b);
        let r = self.fw.append_base(b);
        self.rc.prepend_base(cb);
        r
    }

    #[inline]
    pub fn prepend_base_u8(&mut self, c: u8) -> Base {
        let b = encode_binary_u8(c);
        let cb = complement_base(b);
        let r = self.fw.prepend_base(b);
        self.rc.append_base(cb);
        r
    }

    #[inline]
    pub fn append_base(&mut self, b: Base) -> Base {
        let r = self.fw.append_base(b);
        self.rc.prepend_base(complement_base(b));
        r
    }

    #[inline]
    pub fn prepend_base(&mut self, b: Base) -> Base {
        let r = self.fw.prepend_base(b);
        self.rc.append_base(complement_base(b));
        r
    }

    #[inline]
    pub fn get_canonical_kmer(&self) -> Kmer {
        if self.fw.data < self.rc.data {
            self.fw.clone()
        } else {
            self.rc.clone()
        }
    }

    #[inline]
    pub fn get_canonical_word(&self) -> u64 {
        if self.fw.data < self.rc.data {
            self.fw.data
        } else {
            self.rc.data
        }
    }

    #[inline]
    pub fn get_fw_mer(&self) -> Kmer {
        self.fw.clone()
    }

    #[inline]
    pub fn get_rc_mer(&self) -> Kmer {
        self.rc.clone()
    }

    #[inline]
    pub fn get_fw_word(&self) -> u64 {
        self.fw.data
    }

    #[inline]
    pub fn get_rc_word(&self) -> u64 {
        self.rc.data
    }

    #[inline]
    pub fn get_kmer_equivalency(&self, other: &Kmer) -> MatchType {
        if self.get_fw_word() == other.data {
            MatchType::IdentityMatch
        } else if self.get_rc_word() == other.data {
            MatchType::TwinMatch
        } else {
            MatchType::NoMatch
        }
    }

    #[inline]
    pub fn get_word_equivalency(&self, other: u64) -> MatchType {
        if self.get_fw_word() == other {
            MatchType::IdentityMatch
        } else if self.get_rc_word() == other {
            MatchType::TwinMatch
        } else {
            MatchType::NoMatch
        }
    }
}

impl From<Kmer> for CanonicalKmer {
    #[inline]
    fn from(km: Kmer) -> Self {
        Self {
            rc: km.to_reverse_complement(),
            fw: km,
        }
    }
}

impl From<String> for CanonicalKmer {
    fn from(s: String) -> Self {
        let fk = Kmer::from(s);
        let rk = fk.to_reverse_complement();
        Self { fw: fk, rc: rk }
    }
}

impl From<&str> for CanonicalKmer {
    fn from(s: &str) -> Self {
        Self::from(String::from(s))
    }
}

impl From<&[u8]> for CanonicalKmer {
    fn from(s: &[u8]) -> Self {
        let fw = Kmer::from(s);
        Self {
            fw: fw.clone(),
            rc: fw.to_reverse_complement(),
        }
    }
}

impl From<CanonicalKmer> for String {
    fn from(kmer: CanonicalKmer) -> Self {
        kmer.get_canonical_kmer().into()
    }
}

impl std::fmt::Display for CanonicalKmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let msg: String = self.clone().into();
        write!(f, "{msg}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const K: u8 = 31;

    #[quickcheck]
    fn swap_identity(word: u64) -> bool {
        let mut a = CanonicalKmer::from_u64(word, K);
        let b = a.clone();
        a.swap();
        a.swap();
        a == b
    }

    #[quickcheck]
    fn equivalency(word: u64) -> bool {
        let canon_km = CanonicalKmer::from_u64(word, K);
        let mut canon_km2 = CanonicalKmer::from(canon_km.get_rc_mer());

        let e = canon_km.get_kmer_equivalency(&canon_km2.get_fw_mer());
        let mut res = e == MatchType::TwinMatch;

        canon_km2.swap();
        let e = canon_km.get_kmer_equivalency(&canon_km2.get_fw_mer());
        res &= e == MatchType::IdentityMatch;

        canon_km2.append_base_u8(b'c');
        let e = canon_km.get_kmer_equivalency(&canon_km2.get_fw_mer());
        res &= e == MatchType::NoMatch;
        res
    }

    #[test]
    fn test_from_u64() {
        let km = Kmer::from("acttg");
        let canon_km = CanonicalKmer::from_u64(km.into_u64(), km.len() as u8);

        assert_eq!(canon_km.fw.to_string(), "acttg");
        assert_eq!(canon_km.rc.to_string(), "caagt");
    }

    #[test]
    fn test_from_kmer() {
        let km = Kmer::from("acttg");
        let canon_km = CanonicalKmer::from(km);

        assert_eq!(canon_km.fw.to_string(), "acttg");
        assert_eq!(canon_km.rc.to_string(), "caagt");
    }

    #[test]
    fn test_swap() {
        let mut canon_km = CanonicalKmer::from("acttg");
        assert_eq!(canon_km.fw.to_string(), "acttg");
        assert_eq!(canon_km.rc.to_string(), "caagt");
        canon_km.swap();
        assert_eq!(canon_km.rc.to_string(), "acttg");
        assert_eq!(canon_km.fw.to_string(), "caagt");
    }

    #[test]
    fn test_shift() {
        let mut canon_km = CanonicalKmer::from("acttg");
        canon_km.append_base_u8(b'a');
        assert_eq!(canon_km.fw.to_string(), "cttga");
        assert_eq!(canon_km.rc.to_string(), "tcaag");
        canon_km.prepend_base_u8(b'c');
        assert_eq!(canon_km.rc.to_string(), "caagg");
        assert_eq!(canon_km.fw.to_string(), "ccttg");
    }

    #[test]
    fn test_equivalency() {
        let canon_km = CanonicalKmer::from("acttg");
        let mut canon_km2 = CanonicalKmer::from("caagt");

        let e = canon_km.get_kmer_equivalency(&canon_km2.get_fw_mer());
        assert_eq!(e, MatchType::TwinMatch);

        canon_km2.swap();
        let e = canon_km.get_kmer_equivalency(&canon_km2.get_fw_mer());
        assert_eq!(e, MatchType::IdentityMatch);

        canon_km2.append_base_u8(b'c');
        let e = canon_km.get_kmer_equivalency(&canon_km2.get_fw_mer());
        assert_eq!(e, MatchType::NoMatch);
    }
}
