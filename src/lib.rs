/* mod declaration */
pub mod encoding;
pub mod kmer;
pub mod naive_impl;
pub mod utils;

#[cfg(test)]
extern crate quickcheck;

#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;
