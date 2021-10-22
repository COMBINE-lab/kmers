/* project use */

/* mod declaration */
pub mod naive;
pub mod xor10;

/* public use */
pub use naive::Naive;
pub use xor10::Xor10;

pub trait Encoding<P, const B: usize> {
    fn encode(&self, seq: &[u8]) -> [P; B];

    fn decode(&self, array: [P; B]) -> Vec<u8>;

    fn rev_comp(&self, array: [P; B]) -> [P; B];
}
