/* project use */

/* mod declaration */
pub mod naive;

/* public use */
pub use naive::Naive;

pub trait Encoding<P, const B: usize> {
    fn encode(&self, seq: &[u8]) -> [P; B];

    fn decode(&self, array: [P; B]) -> Vec<u8>;

    fn rev_comp(&self, array: [P; B]) -> [P; B];
}
