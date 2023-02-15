use std::collections::VecDeque;
use std::hash::BuildHasher;

use super::super::hash::hash_one;
use super::*;

#[derive(Clone, Debug)]
struct DQMer {
    pub lmer: u64,
    pub pos: usize,
    pub hash: u64,
}

impl DQMer {
    fn new(lmer: u64, pos: usize, hash: u64) -> Self {
        Self { lmer, pos, hash }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MappedMinimizer {
    word: u64,  // u64 representation
    pub pos: usize, // position in sequence
}

impl MappedMinimizer {
    pub fn as_u64(&self) -> u64{
        self.word
    }
}

impl MappedMinimizer {
    pub fn new(lmer: u64, pos: usize) -> Self {
        Self { word: lmer, pos }
    }
}

pub struct SeqVecMinimizerIter<'a, T: BuildHasher> {
    dq: VecDeque<DQMer>,
    k: usize,
    w: usize, // or "L"
    curr_km_i: usize,
    sv: SeqVectorSlice<'a>,
    hash_seed: T,
}

impl<'a, T: BuildHasher> SeqVecMinimizerIter<'a, T> {
    // Dequeue invariant:
    // Q = [ (Li, pi) ... (Lj, pj) ] at iteration ii
    //   front                    back

    // - each L_i is an L-mer in the (k-1)-prefix of kmer_ii
    // - where each L_i can be a candidate minimizer in any suffix of kmer_ii
    //   and thus could be a minimizer in any kmer succeeding kmer_ii with at least a L-mer overlap
    // - i.e.:
    //    - h(L_i) < h(L_j) iff i < j
    //    - pj - pi
    //

    #[inline]
    fn enqueue_dqmer(&mut self, dqmer: DQMer) {
        // Could be a while loop to be safe...
        // but we should never have more than one minimizer fall out of the window
        // since only one successive lmer is added to the queue

        if let Some(frontmer) = self.dq.front() {
            if frontmer.pos < self.curr_km_i {
                self.dq.pop_front();
            }
        }

        while let Some(backmer) = self.dq.back() {
            if backmer.hash <= dqmer.hash {
                break;
            } else {
                self.dq.pop_back();
            }
        }

        self.dq.push_back(dqmer);
    }

    #[inline]
    fn next_dqmer(&mut self) -> DQMer {
        // return last dqmer of curr_km_ii-th kmer
        let pos = self.curr_km_i + self.k - self.w;
        let lmer = self.sv.get_kmer_u64(pos, self.w);
        let hash = hash_one(&self.hash_seed, lmer);
        DQMer::new(lmer, pos, hash)
    }

    #[inline]
    fn n_kmers(&self) -> usize {
        self.sv.len() - self.k + 1
    }

    pub fn new(sv: SeqVectorSlice<'a>, k: usize, w: usize, hash_seed: T) -> Self {
        // Insert lmers of the k-1 prefix
        assert!(sv.len() >= k);
        let dq = VecDeque::with_capacity(k - w + 1);

        let mut iter = Self {
            dq,
            k,
            w,
            hash_seed,
            sv: sv.clone(),
            curr_km_i: 0,
        };

        for i in 0..(k - w) {
            let lmer = sv.get_kmer_u64(i, w);
            let hash = hash_one(&iter.hash_seed, lmer);

            let dqmer = DQMer { lmer, pos: i, hash };

            iter.enqueue_dqmer(dqmer)
        }

        iter
    }
}

impl<T: BuildHasher> Iterator for SeqVecMinimizerIter<'_, T> {
    type Item = MappedMinimizer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_km_i < self.n_kmers() {
            let dqmer = self.next_dqmer();
            self.enqueue_dqmer(dqmer);
            let dqmer = self.dq.front().unwrap();
            let mmer = MappedMinimizer {
                word: dqmer.lmer,
                pos: dqmer.pos,
            };
            self.curr_km_i += 1;
            Some(mmer)
        } else {
            None
        }
    }
}

// Minimizer iterator does not know how many minimizers there are or its length
// impl<T: BuildHasher> ExactSizeIterator for SeqVecMinimizerIter<'_, T> {
//     fn len(&self) -> usize {
//         self.sv.len() - self.curr_km_i - self.k + 1
//     }
// }

#[cfg(test)]
mod test {
    use std::collections::hash_map::RandomState;
    use std::collections::VecDeque;

    use crate::naive_impl::hash::LexHasherState;

    use super::*;

    fn dqmers_from_hashes(hashes: &[u64]) -> Vec<DQMer> {
        hashes
            .iter()
            .enumerate()
            .map(|(i, h)| DQMer::new(0, i, *h))
            .collect()
    }

    fn iter_dq_hashes<T: BuildHasher>(iter: &SeqVecMinimizerIter<T>) -> Vec<u64> {
        iter.dq.iter().map(|x| x.hash).collect()
    }

    #[test]
    fn enqueue_dqmer() {
        let sv = SeqVector::from(b"");
        let sv = sv.as_slice();

        let dq = VecDeque::new();

        let (k, w) = (4, 2);
        let hashes = vec![2, 1, 0, 0, 3, 4, 2];

        let mut iter = SeqVecMinimizerIter {
            sv,
            dq,
            k,
            w,
            curr_km_i: 0,
            hash_seed: RandomState::new(),
        };

        let dqmers = dqmers_from_hashes(&hashes);

        iter.enqueue_dqmer(dqmers[0].clone());
        assert_eq!(iter_dq_hashes(&iter), vec![2]);

        iter.enqueue_dqmer(dqmers[1].clone());
        assert_eq!(iter_dq_hashes(&iter), vec![1]);

        iter.enqueue_dqmer(dqmers[2].clone());
        assert_eq!(iter_dq_hashes(&iter), vec![0]);
        iter.curr_km_i = 1;

        iter.enqueue_dqmer(dqmers[3].clone());
        assert_eq!(iter_dq_hashes(&iter), vec![0, 0]);
        iter.curr_km_i = 2;

        iter.enqueue_dqmer(dqmers[4].clone());
        assert_eq!(iter_dq_hashes(&iter), vec![0, 0, 3]);
        iter.curr_km_i = 3;

        iter.enqueue_dqmer(dqmers[5].clone());
        assert_eq!(iter_dq_hashes(&iter), vec![0, 3, 4]);
        iter.curr_km_i = 4;

        iter.enqueue_dqmer(dqmers[6].clone());
        assert_eq!(iter_dq_hashes(&iter), vec![2]);
        iter.curr_km_i = 5;
    }

    #[test]
    fn leftmost_mmer() {
        let sv = SeqVector::from(b"AAAAAAA");
        let iter = SeqVecMinimizerIter::new(sv.as_slice(), 5, 3, RandomState::new());

        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::new(0, 0),
                MappedMinimizer::new(0, 1),
                MappedMinimizer::new(0, 2),
            ]
        )
    }

    #[test]
    fn mmers0() {
        let sv = SeqVector::from(b"AAACAAA");
        let iter = SeqVecMinimizerIter::new(sv.as_slice(), 6, 3, LexHasherState::new(6));

        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![MappedMinimizer::new(0, 0), MappedMinimizer::new(0, 4),]
        )
    }

    #[test]
    fn mmers1() {
        let sv = SeqVector::from(b"AACCAAA");
        let iter = SeqVecMinimizerIter::new(sv.as_slice(), 5, 3, LexHasherState::new(5));

        let mmers: Vec<MappedMinimizer> = iter.collect();

        let aac = 0b010000;
        let acc = 0b010100;
        let aaa = 0b000000;
        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::new(aac, 0),
                MappedMinimizer::new(acc, 1),
                MappedMinimizer::new(aaa, 4),
            ]
        )
    }

    #[test]
    fn mmers2() {
        let sv = SeqVector::from(b"CACACACCAC");
        // let bh = RandomState::new();
        let bh = LexHasherState::new(3);
        let iter = SeqVecMinimizerIter::new(sv.as_slice(), 7, 3, bh);

        let mmers: Vec<MappedMinimizer> = iter.collect();

        // let aac = 0b010000;
        let aca = 0b000100;
        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::new(aca, 1),
                MappedMinimizer::new(aca, 1),
                MappedMinimizer::new(aca, 3),
                MappedMinimizer::new(aca, 3),
            ]
        )
    }
}
