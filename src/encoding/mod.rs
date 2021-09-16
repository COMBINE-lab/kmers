/* project use */

/* mod declaration */
pub mod naive;

pub trait Encoder<P, const B: usize> {
    fn encode(&self, seq: &[u8]) -> [P; B];
}
