use std::mem;

struct Kmer<const K: u16, const B: usize> {
    array: [u64; B],
}

impl<const K: u16, const B: usize> Kmer<K, B> {
    /// construct a new empty k-mer 
    /// it is constructed as poly-A by default
    pub fn new() -> Self {
        Kmer {
            array: [0;B]
        }
    }

    /// returns the number of bytes used for the 
    /// storage of this k-mer
    pub fn num_bytes(&self) -> usize {
        std::mem::size_of::<u64>() * self.array.len()
    }

    /// returns the value of k for this k-mer
    pub fn k(&self) -> usize {
        K as usize
    }
}

#[cfg(test)]
mod tests {
    use crate::Kmer;
    #[test]
    fn it_works() {
        let km = Kmer::<15,1>::new();
        assert_eq!(km.num_bytes(), 8);
        assert_eq!(km.k(), 15);
    }
}
