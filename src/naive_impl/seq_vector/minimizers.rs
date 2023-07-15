use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::marker::PhantomData;

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct MappedMinimizer {
    word: u64,  // u64 representation
    pos: usize, // position in sequence
}

impl MappedMinimizer {
    pub fn as_u64(&self) -> u64 {
        self.word
    }

    pub fn pos(&self) -> usize {
        self.pos
    }

    pub fn from_parts(lmer: u64, pos: usize) -> Self {
        Self { word: lmer, pos }
    }

    pub fn from_seq<T>(seq: T, pos: usize) -> Self
    where
        crate::naive_impl::Kmer: From<T>,
    {
        let word = crate::naive_impl::Kmer::from(seq).into_u64();
        Self { word, pos }
    }
}

#[derive(PartialEq, Debug)]
pub struct LeftMin;

#[derive(PartialEq, Debug)]
pub struct RightMin;

#[derive(PartialEq, Debug)]
pub struct HashedMinimizer<T> {
    pub word: u64,
    pub hash: u64,
    pub pos: usize,
    pub phantom: PhantomData<T>,
}

impl<T> HashedMinimizer<T> {
    fn new(word: u64, hash: u64, pos: usize) -> Self {
        Self {
            word,
            hash,
            pos,
            phantom: PhantomData::default(),
        }
    }

    fn to_mapped_minimizer(&self) -> MappedMinimizer {
        MappedMinimizer {
            word: self.word,
            pos: self.pos,
        }
    }
}

impl PartialOrd for HashedMinimizer<LeftMin> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        let ord = match self.hash.cmp(&rhs.hash) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.pos.cmp(&rhs.pos),
        };
        Some(ord)
    }
}

impl PartialOrd for HashedMinimizer<RightMin> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        let ord = match self.hash.cmp(&rhs.hash) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.pos.cmp(&rhs.pos).reverse(),
        };
        Some(ord)
    }
}

use super::super::hash::hash_one;
use super::SeqVectorSlice;
use std::hash::BuildHasher;

struct HashedMinimizerQueue<OrdT> {
    q: VecDeque<HashedMinimizer<OrdT>>,
}

impl<OrdT> HashedMinimizerQueue<OrdT>
where
    HashedMinimizer<OrdT>: PartialOrd,
{
    pub fn with_capacity(capacity: usize) -> Self {
        let q = VecDeque::with_capacity(capacity);

        Self { q }
    }
    pub fn front(&self) -> Option<&HashedMinimizer<OrdT>> {
        self.q.front()
    }

    // Queue maintains the minimizers of all suffixes of a k-mer window
    pub fn eat(&mut self, mmer: HashedMinimizer<OrdT>, min_pos: usize) {
        if let Some(frontmer) = self.q.front() {
            // remove wmer that falls out of window
            if frontmer.pos < min_pos {
                self.q.pop_front();
            }
        }

        while let Some(backmer) = self.q.back() {
            // update suffix minimizers
            if backmer <= &mmer {
                break;
            } else {
                self.q.pop_back();
            }
        }

        self.q.push_back(mmer);
    }
}

pub type MinimizerIterLeftMin<'a, T> = MinimizerIter<'a, T, LeftMin>;
pub type MinimizerIterRightMin<'a, T> = MinimizerIter<'a, T, RightMin>;

pub struct MinimizerIter<'a, T: BuildHasher, OrdT> {
    q: HashedMinimizerQueue<OrdT>,
    k: usize,
    w: usize, // or "L"
    curr_km_i: usize,
    sv: SeqVectorSlice<'a>,
    hash_seed: T,
}

impl<'a, T: BuildHasher, OrdT> MinimizerIter<'a, T, OrdT>
where
    HashedMinimizer<OrdT>: PartialOrd,
{
    #[inline]
    fn get_mmer(&self, pos: usize) -> HashedMinimizer<OrdT> {
        let lmer = self.sv.get_kmer_u64(pos, self.w);
        let hash = hash_one(&self.hash_seed, lmer);
        HashedMinimizer::new(lmer, hash, pos)
    }

    #[inline]
    fn next_mmer(&self) -> HashedMinimizer<OrdT> {
        // return last HashedMinimizer of curr_km_ii-th kmer
        let pos = self.curr_km_i + self.k - self.w;
        self.get_mmer(pos)
    }

    #[inline]
    fn n_kmers(&self) -> usize {
        self.sv.len() - self.k + 1
    }

    pub fn new(sv: SeqVectorSlice<'a>, k: usize, w: usize, hash_seed: T) -> Self {
        // Insert lmers of the k-1 prefix
        assert!(sv.len() >= k);

        let q = HashedMinimizerQueue::with_capacity(k - w + 1);
        let mut iter = Self {
            q,
            k,
            w,
            hash_seed,
            sv: sv.clone(),
            curr_km_i: 0,
        };

        for i in 0..(k - w) {
            let mmer = iter.get_mmer(i);
            iter.q.eat(mmer, 0)
        }

        iter
    }

    #[inline]
    pub fn eat_next_mmer(&mut self) {
        let mmer = self.next_mmer();
        self.q.eat(mmer, self.curr_km_i);
    }
}

impl<T: BuildHasher, OrdT> Iterator for MinimizerIter<'_, T, OrdT>
where
    HashedMinimizer<OrdT>: PartialOrd,
{
    type Item = MappedMinimizer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_km_i < self.n_kmers() {
            self.eat_next_mmer();
            let mmer = self.q.front().unwrap();
            let mmer = MappedMinimizer {
                word: mmer.word,
                pos: mmer.pos,
            };
            self.curr_km_i += 1;
            Some(mmer)
        } else {
            None
        }
    }
}

pub struct CanonicalMinimizerIter<'a, T: BuildHasher> {
    fwq: HashedMinimizerQueue<LeftMin>,
    rcq: HashedMinimizerQueue<RightMin>,

    k: usize,
    w: usize, // or "L"
    curr_km_i: usize,
    sv: SeqVectorSlice<'a>,
    hash_seed: T,
}

type MinimizerPair = (HashedMinimizer<LeftMin>, HashedMinimizer<RightMin>);
impl<'a, T: BuildHasher> CanonicalMinimizerIter<'a, T> {
    pub fn new(sv: SeqVectorSlice<'a>, k: usize, w: usize, hash_seed: T) -> Self {
        // Insert lmers of the k-1 prefix
        assert!(sv.len() >= k);

        let fwq = HashedMinimizerQueue::with_capacity(k - w + 1);
        let rcq = HashedMinimizerQueue::with_capacity(k - w + 1);

        let mut iter = Self {
            fwq,
            rcq,
            k,
            w,
            hash_seed,
            sv: sv.clone(),
            curr_km_i: 0,
        };

        for i in 0..(k - w) {
            let (fw_mmer, rc_mmer) = iter.get_mmer_pair(i);
            iter.fwq.eat(fw_mmer, 0);
            iter.rcq.eat(rc_mmer, 0);
        }

        iter
    }

    #[inline]
    fn get_mmer_pair(&self, pos: usize) -> MinimizerPair {
        let fw_mmer = self.sv.get_kmer(pos, self.w);
        let rc_mmer = fw_mmer.to_reverse_complement().into_u64();
        let fw_mmer = fw_mmer.into_u64();

        let fw_hash = hash_one(&self.hash_seed, fw_mmer);
        let rc_hash = hash_one(&self.hash_seed, rc_mmer);

        let fw_mmer = HashedMinimizer::new(fw_mmer, fw_hash, pos);
        let rc_mmer = HashedMinimizer::new(rc_mmer, rc_hash, pos);
        (fw_mmer, rc_mmer)
    }

    #[inline]
    fn next_mmer_pair(&self) -> MinimizerPair {
        // return last HashedMinimizer of curr_km_ii-th kmer
        let pos = self.curr_km_i + self.k - self.w;
        self.get_mmer_pair(pos)
    }

    #[inline]
    fn n_kmers(&self) -> usize {
        self.sv.len() - self.k + 1
    }

    #[inline]
    pub fn eat_next_mmer_pair(&mut self) {
        let (fw_mmer, rc_mmer) = self.next_mmer_pair();

        self.fwq.eat(fw_mmer, self.curr_km_i);
        self.rcq.eat(rc_mmer, self.curr_km_i);
    }
}

impl<T: BuildHasher> Iterator for CanonicalMinimizerIter<'_, T> {
    type Item = MappedMinimizer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_km_i < self.n_kmers() {
            self.eat_next_mmer_pair(); // add next minimizer candidates in fw and rc

            let km = self.sv.get_kmer(self.curr_km_i, self.k);

            let mmer = if km.is_canonical() {
                self.fwq.front().unwrap().to_mapped_minimizer()
            } else {
                self.rcq.front().unwrap().to_mapped_minimizer()
            };

            self.curr_km_i += 1;
            Some(mmer)
        } else {
            None
        }
    }
}

pub struct CanonicalSuperKmerIterator<'a, T: BuildHasher> {
    minimizers: CanonicalMinimizerIter<'a, T>,
    k: usize,
    w: usize, // or "L"
    curr_km_i: usize,
    sv: SeqVectorSlice<'a>,
    next_mmer: Option<MappedMinimizer>, // next mmer to consume.
}

#[derive(Clone, PartialEq, Debug)]
pub struct CanonicalSuperKmerOcc {
    mmer: MappedMinimizer,
    start: usize,
    n_kmers: u8, // a little optimization to use u8 insead of usize
                 // for external memory construction
}

impl CanonicalSuperKmerOcc {
    pub fn from_parts(mmer: MappedMinimizer, start: usize, n_kmers: u8) -> Self {
        Self {
            mmer,
            start,
            n_kmers,
        }
    }

    pub fn inc_pos(&mut self, offset: usize) {
        // increase the position of the super-k-mer by given offset.
        self.mmer.pos += offset;
        self.start += offset;
    }

    pub fn dec_pos(&mut self, offset: usize) {
        // decrease the position of the super-k-mer by a given offset.
        assert!(offset <= self.start);
        self.mmer.pos -= offset;
        self.start -= offset;
    }

    #[inline]
    pub fn mmer_word(&self) -> u64 {
        self.mmer.word
    }

    #[inline]
    pub fn mmer_pos(&self) -> usize {
        self.mmer.pos
    }

    #[inline]
    pub fn mmer_offset(&self) -> usize {
        self.mmer_pos() - self.start_pos()
    }

    #[inline]
    pub fn start_pos(&self) -> usize {
        self.start
    }

    #[inline]
    pub fn n_kmers(&self) -> usize {
        self.n_kmers as usize
    }
}

impl<'a, T: BuildHasher> CanonicalSuperKmerIterator<'a, T> {
    pub fn new(sv: SeqVectorSlice<'a>, k: usize, w: usize, hash_seed: T) -> Self {
        let mut minimizers = CanonicalMinimizerIter::new(sv.clone(), k, w, hash_seed);
        let next_mmer = minimizers.next();

        Self {
            minimizers: minimizers,
            sv,
            k,
            w,
            curr_km_i: 0,
            next_mmer,
        }
    }
}

impl<'a, T: BuildHasher> Iterator for CanonicalSuperKmerIterator<'a, T> {
    type Item = CanonicalSuperKmerOcc;

    fn next(&mut self) -> Option<Self::Item> {
        let curr_mmer = self.next_mmer.as_ref()?.clone();
        let start_pos = self.curr_km_i;
        let mut n_kmers = 1;

        loop {
            // update the current minimiizer
            self.next_mmer = self.minimizers.next();
            self.curr_km_i += 1;

            let keep_searching = self.next_mmer.as_ref().is_some_and(|next| {
                (next.pos == curr_mmer.pos) && (next.word == curr_mmer.word)
            });

            if keep_searching {
                n_kmers += 1;
            } else {
                // Either curr_mmer is diff from next, or curr_mmer is none
                // and we are at last super-kmer.
                return Some(CanonicalSuperKmerOcc {
                    mmer: curr_mmer,
                    start: start_pos,
                    n_kmers,
                });
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::SeqVector;
    use super::*;
    use crate::naive_impl::hash::LexHasherState;
    use std::collections::hash_map::RandomState;

    #[test]
    fn leftmost_mmer() {
        let sv = SeqVector::from(b"AAAAAAA");
        let iter = MinimizerIterLeftMin::new(sv.as_slice(), 5, 3, RandomState::new());

        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"AAA", 0),
                MappedMinimizer::from_seq(b"AAA", 1),
                MappedMinimizer::from_seq(b"AAA", 2),
            ]
        )
    }

    #[test]
    fn rightmost_mmer() {
        let sv = SeqVector::from(b"AAAAAAA");
        let iter = MinimizerIterRightMin::new(sv.as_slice(), 5, 3, RandomState::new());

        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"AAA", 2),
                MappedMinimizer::from_seq(b"AAA", 3),
                MappedMinimizer::from_seq(b"AAA", 4),
            ]
        )
    }

    #[test]
    fn mmers0() {
        let sv = SeqVector::from(b"AAACAAA");
        let iter = MinimizerIterLeftMin::new(sv.as_slice(), 6, 3, LexHasherState::new(6));

        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"AAA", 0),
                MappedMinimizer::from_seq(b"AAA", 4),
            ]
        )
    }

    #[test]
    fn mmers1() {
        let sv = SeqVector::from(b"AACCAAA");
        let iter = MinimizerIterLeftMin::new(sv.as_slice(), 5, 3, LexHasherState::new(5));
        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"AAC", 0),
                MappedMinimizer::from_seq(b"ACC", 1),
                MappedMinimizer::from_seq(b"AAA", 4),
            ]
        )
    }

    #[test]
    fn mmers2() {
        let sv = SeqVector::from(b"CACACACCAC");
        // let bh = RandomState::new();
        let bh = LexHasherState::new(3);
        let iter = MinimizerIterLeftMin::new(sv.as_slice(), 7, 3, bh);

        let mmers: Vec<MappedMinimizer> = iter.collect();

        // let aac = 0b010000;
        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"ACA", 1),
                MappedMinimizer::from_seq(b"ACA", 1),
                MappedMinimizer::from_seq(b"ACA", 3),
                MappedMinimizer::from_seq(b"ACA", 3),
            ]
        )
    }
}

#[cfg(test)]
mod test_canonical {
    use super::super::SeqVector;
    use super::*;
    use crate::naive_impl::hash::LexHasherState;
    use std::collections::hash_map::RandomState;

    // Note - Defining minimizers and thier positions on canonical k-mers
    // - mini(g*) = mini( min(g, g') )
    //   breaking ties by taking leftmost (least index) w-mer on the canonical sequence g*.
    // - if g != min(g, g'), i.e. is not canonical, the position of the minimizer
    //   occurrence is the _rightmost_ occurrence of mini(g')' = x'.

    #[test]
    fn leftmost_on_canon() {
        let sv = SeqVector::from(b"AAAAAAA");
        let iter = CanonicalMinimizerIter::new(sv.as_slice(), 5, 3, RandomState::new());

        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"AAA", 0),
                MappedMinimizer::from_seq(b"AAA", 1),
                MappedMinimizer::from_seq(b"AAA", 2),
            ]
        )
    }

    #[test]
    fn rightmost_on_rc() {
        let sv = SeqVector::from(b"TTTTTTT");
        let iter = CanonicalMinimizerIter::new(sv.as_slice(), 5, 3, RandomState::new());

        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"AAA", 2),
                MappedMinimizer::from_seq(b"AAA", 3),
                MappedMinimizer::from_seq(b"AAA", 4),
            ]
        )
    }

    #[test]
    fn break_ties_on_canonical() {
        // Check that ties are broken properly on canonical sequence
        // minimizer is the leftmost minimum w-mer on canonical sequence
        // Given the seq TAAATTTC
        // 1. The first kmer  TAAA[TTT] canonicalized is [aaa]ttta with  minmizer in [...] @ idx == 4
        // 2. The second kmer  [AAA]TTTC is already canonical with minimizer AAA @ idx 1

        let (k, w) = (7, 3);
        let sv = SeqVector::from(b"TAAATTTC");
        let iter = CanonicalMinimizerIter::new(sv.as_slice(), k, w, LexHasherState::new(w));
        let mmers: Vec<MappedMinimizer> = iter.collect();

        assert_eq!(
            mmers,
            vec![
                MappedMinimizer::from_seq(b"AAA", 4),
                MappedMinimizer::from_seq(b"AAA", 1),
            ]
        )
    }

    #[test]
    fn super_kmers() {
        // Super k-mers example
        // K=7, t = 3

        //                    : A G G G A A A G A A
        //    AGGGAAA(tttccct): A G G G[A A A]       // s0
        //    GGGAAAG(ctttccc):  [G G G]A A A G      // s1
        //    GGAAAGA(tctttgg):     G G[A A A]G A    // s2
        //    GAAAGAA(ttctttc):       G[A A A]G A A  // s2

        let (k, w) = (7, 3);
        let sv = SeqVector::from(b"AGGGAAAGAA");
        let iter = CanonicalSuperKmerIterator::new(sv.as_slice(), k, w, LexHasherState::new(w));
        let skms: Vec<CanonicalSuperKmerOcc> = iter.collect();

        assert_eq!(skms.len(), 3);

        assert_eq!(
            skms[0],
            CanonicalSuperKmerOcc {
                mmer: MappedMinimizer::from_seq(b"AAA", 4),
                start: 0,
                n_kmers: 1,
            }
        );

        assert_eq!(
            skms[1],
            CanonicalSuperKmerOcc {
                mmer: MappedMinimizer::from_seq(b"CCC", 1),
                start: 1,
                n_kmers: 1,
            }
        );

        assert_eq!(
            skms[2],
            CanonicalSuperKmerOcc {
                mmer: MappedMinimizer::from_seq(b"AAA", 4),
                start: 2,
                n_kmers: 2,
            }
        )
    }

    #[test]
    fn super_kmers2() {
        // Super k-mers example with minimizers in same position but as rc and different
        // minimizer seqs on adjacent kmers.
        // k=5,t=3
        // G G C T T A
        // G G[C T T] (aagcc)  => with  mini seq aag since rc is canonical
        //   G[C T T]A(taagc)  => with mini seq CTT since fw is canonical

        let (k, w) = (5, 3);
        let sv = SeqVector::from(b"GGCTTA");
        let iter = CanonicalSuperKmerIterator::new(sv.as_slice(), k, w, LexHasherState::new(w));
        let skms: Vec<CanonicalSuperKmerOcc> = iter.collect();

        assert_eq!(skms.len(), 2);

        assert_eq!(
            skms[0],
            CanonicalSuperKmerOcc {
                mmer: MappedMinimizer::from_seq(b"aag", 2),
                start: 0,
                n_kmers: 1,
            }
        );

        assert_eq!(
            skms[1],
            CanonicalSuperKmerOcc {
                mmer: MappedMinimizer::from_seq(b"CTT", 2),
                start: 1,
                n_kmers: 1,
            }
        );
    }
}
