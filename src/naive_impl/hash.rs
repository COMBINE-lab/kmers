use super::Kmer;
use std::hash::{BuildHasher, Hash, Hasher};

impl Hash for Kmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.data);
    }
}

pub fn hash_one<H, T>(state: &H, x: T) -> u64
where
    H: BuildHasher,
    T: Hash,
    H: Sized,
    H::Hasher: Hasher,
{
    let mut hasher = state.build_hasher();
    x.hash(&mut hasher);
    hasher.finish()
}

#[derive(Clone, Debug, PartialEq)]
pub struct LexHasherState(usize);

impl LexHasherState {
    pub fn new(k: usize) -> Self {
        Self(k)
    }
}

impl BuildHasher for LexHasherState {
    type Hasher = LexHasher;
    fn build_hasher(&self) -> Self::Hasher {
        LexHasher::new(self.0)
    }
}

// TODO impl with debug assertions with #[cfg(debug_assertions)]
// like https://github.com/paritytech/nohash-hasher/.../lib.rs#L106
pub struct LexHasher {
    state: u64,
    k: usize,
}

impl LexHasher {
    pub fn new(k: usize) -> Self {
        Self { k, state: 0 }
    }
}

impl Hasher for LexHasher {
    fn write(&mut self, _: &[u8]) {
        unimplemented!("Hash with write_u64");
    }

    fn finish(&self) -> u64 {
        self.state
    }

    fn write_u64(&mut self, word: u64) {
        // reverse and shift
        let mut res = word;
        res = (res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2;
        res = (res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4;
        res = (res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8;
        res = (res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16;
        res = (res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32;

        res >>= (32 - self.k) * 2;
        self.state = res;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn new() {
        let s = LexHasher::new(3);
        assert_eq!(s.k, 3);
    }
    #[test]
    fn lex_order() {
        let seed = LexHasherState::new(3);

        let aaa = Kmer::from(b"aaa");
        let h1 = hash_one(&seed, aaa);
        assert_eq!(h1, 0);

        let aac = Kmer::from(b"aac");

        let h2 = hash_one(&seed, aac);
        assert!(h1 < h2);

        assert_eq!(h2, 0b00001);

        let cac = hash_one(&seed, Kmer::from(b"cac"));
        let caa = hash_one(&seed, Kmer::from(b"caa"));

        assert!(caa < cac);
        assert_eq!(caa, 0b010000);
        assert_eq!(cac, 0b010001);
    }
}
