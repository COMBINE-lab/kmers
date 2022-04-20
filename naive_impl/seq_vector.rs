use serde::{Deserialize, Serialize};
use simple_sds::int_vector::IntVector;
use simple_sds::ops::Vector;
use simple_sds::raw_vector::{AccessRaw, RawVector};

use super::Kmer;
use crate::serde_ext;

#[derive(Clone, Debug, Serialize, Deserialize, Eq, PartialEq)]
pub struct SeqVector {
    // k: u8,
    #[serde(with = "serde_ext")]
    data: RawVector,
}

impl SeqVector {
    pub fn len(&self) -> usize {
        self.data.len() / 2
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn get_kmer(&self, pos: usize, k: u8) -> Kmer {
        Kmer::from_u64(self.get_kmer_u64(pos, k), k)
    }

    pub fn get_kmer_u64(&self, pos: usize, k: u8) -> u64 {
        assert!(pos < self.len());
        unsafe { self.data.int(pos * 2, k as usize * 2) }
    }

    pub fn get_base(&self, pos: usize) -> u64 {
        self.get_kmer_u64(pos, 1)
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
