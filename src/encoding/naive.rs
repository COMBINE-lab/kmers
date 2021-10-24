//! Naive implementation of encoding
//! Not the best implementation but you can use any encoding
//!
//! This implementation use encoding A -> 00, C -> 01, T -> 10 and G -> 11 to perform conversion to any other encoding.
//! We use this encoding because it's easy to convert ASCII value of nucleotide in this encoding, by perform a bit shift on rigth and keep two lower bits.
//! This work for lower and upper case.

/* crate use */
use bit_field::BitArray as _;

/* Private function */

/// Convert ASCII nucleotide in internal encoding
fn nuc2internal(nuc: u8) -> u8 {
    (nuc as u8 >> 1) & 0b11
}

/// Lookup table usefull to convert internal encoding in ASCII
const INTERNAL2NUC: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Just a wrapper around INTERNAL2NUC lookup table
#[inline]
fn internal2nuc(internal: u8) -> u8 {
    INTERNAL2NUC[internal as usize]
}

/// Function to convert an encoding in a reverse encoding use to perform 2 bits to 8 bits operation
#[inline]
const fn rev_encoding(encoding: u8) -> u8 {
    let mut rev = 0;

    // Sorry for black magic bits operation
    rev ^= 0b00 << (6 - ((encoding >> 6) * 2));
    rev ^= 0b01 << (6 - (((encoding >> 4) & 0b11) * 2));
    rev ^= 0b10 << (6 - (((encoding >> 2) & 0b11) * 2));
    rev ^= 0b11 << (6 - ((encoding & 0b11) * 2));

    rev
}

/* Public interface */

/// Enumeration of all possible encoding:
/// - the first nucleotide is equal to 00
/// - the second nucleotide is equal to 01
/// - the third nucleotide is equal to 10
/// - the last nucleotide is equal to 11
#[derive(Clone, Copy, Debug)]
pub enum Naive {
    ACTG = 0b_00_01_10_11,
    ACGT = 0b_00_01_11_10,
    ATCG = 0b_00_10_01_11,
    ATGC = 0b_00_11_01_10,
    AGCT = 0b_00_10_11_01,
    AGTC = 0b_00_11_10_01,
    CATG = 0b_01_00_10_11,
    CAGT = 0b_01_00_11_10,
    CTAG = 0b_10_00_01_11,
    CTGA = 0b_11_00_01_10,
    CGAT = 0b_10_00_11_01,
    CGTA = 0b_11_00_10_01,
    TACG = 0b_01_10_00_11,
    TAGC = 0b_01_11_00_10,
    TCAG = 0b_10_01_00_11,
    TCGA = 0b_11_01_00_10,
    TGAC = 0b_10_11_00_01,
    TGCA = 0b_11_10_00_01,
    GACT = 0b_01_10_11_00,
    GATC = 0b_01_11_10_00,
    GCAT = 0b_10_01_11_00,
    GCTA = 0b_11_01_10_00,
    GTAC = 0b_10_11_01_00,
    GTCA = 0b_11_10_01_00,
}

impl Naive {
    /// Convert nucleotide in encoding corresponding 2 bits
    pub(crate) fn nuc2bits<P>(&self, nuc: u8) -> P
    where
        P: crate::utils::Data,
    {
        let index = 6 - nuc2internal(nuc) * 2;

        P::from((*self as u8 >> index) & 0b11)
    }

    /// Convert nucleotide encode on 2 bits in 8 bits encoding
    pub(crate) fn bits2nuc<P>(&self, bits: P) -> u8
    where
        P: crate::utils::Data,
    {
        // I forget how this work but it's works
        let rev_encoding = rev_encoding(*self as u8);
        internal2nuc((rev_encoding >> (6 - (bits.to_u8() & 0b11) * 2)) & 0b11)
    }

    /// Get the complement of a nucleotide encode in 2 bits
    pub(crate) fn complement<P>(&self, bits: P) -> P
    where
        P: crate::utils::Data,
    {
        let rev_encoding = rev_encoding(*self as u8);

        let internal = (rev_encoding >> (6 - (bits.to_u8() & 0b11) * 2)) & 0b11;

        let comp_internal = (internal ^ 0b10) & 0b11;

        P::from((*self as u8 >> (6 - comp_internal * 2)) & 0b11)
    }
}

impl<P, const B: usize> super::Encoding<P, B> for Naive
where
    P: crate::utils::Data,
{
    fn encode(&self, seq: &[u8]) -> [P; B] {
        let mut array: [P; B] = unsafe { [std::mem::zeroed(); B] };

        for (idx, nuc) in seq.iter().enumerate() {
            array.set_bits(idx * 2..=idx * 2 + 1, self.nuc2bits(*nuc));
        }

        array
    }

    fn decode(&self, array: [P; B]) -> Vec<u8> {
        let mut seq = Vec::with_capacity(B * P::BIT_LENGTH);

        for idx in 0..array.len() * P::BIT_LENGTH / 2 {
            let value = array.get_bits(idx * 2..=idx * 2 + 1);

            seq.push(self.bits2nuc(value));
        }

        seq
    }

    fn rev_comp(&self, mut array: [P; B]) -> [P; B] {
        let mut i = 0;
        let mut j = array.len() * P::BIT_LENGTH - 2;

        while i < j {
            let comp_i = self.complement(array.get_bits(i..i + 2));
            let comp_j = self.complement(array.get_bits(j..j + 2));

            array.set_bits(i..i + 2, comp_j);
            array.set_bits(j..j + 2, comp_i);

            i += 2;
            j -= 2;
        }

        // No need to manage odd case,
        // array is always a sum of power of 2 bits so its always even

        array
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bit_field::BitField as _;

    use crate::kmer;

    use crate::encoding::Encoding as _;

    #[test]
    fn one_base_all_encoding() {
        macro_rules! one_base2bits {
	    ($($ty:expr), *) => (
		$(
		    let encoder = $ty;

		    assert_eq!(encoder.nuc2bits::<u8>(b'A'), (encoder as u8).get_bits(6..8));
		    assert_eq!(encoder.nuc2bits::<u8>(b'C'), (encoder as u8).get_bits(4..6));
		    assert_eq!(encoder.nuc2bits::<u8>(b'T'), (encoder as u8).get_bits(2..4));
		    assert_eq!(encoder.nuc2bits::<u8>(b'G'), (encoder as u8).get_bits(0..2));
		)*
	    )
	}

        one_base2bits!(
            Naive::ACTG,
            Naive::ACGT,
            Naive::ATCG,
            Naive::ATGC,
            Naive::AGCT,
            Naive::AGTC,
            Naive::CATG,
            Naive::CAGT,
            Naive::CTAG,
            Naive::CTGA,
            Naive::CGAT,
            Naive::CGTA,
            Naive::TACG,
            Naive::TAGC,
            Naive::TCAG,
            Naive::TCGA,
            Naive::TGAC,
            Naive::TGCA,
            Naive::GACT,
            Naive::GATC,
            Naive::GCAT,
            Naive::GCTA,
            Naive::GTAC,
            Naive::GTCA
        );
    }

    #[test]
    fn one_base_all_decoding() {
        macro_rules! bits2one_base {
	    ($($ty:expr), *) => (
		$(
		    let encoder = $ty;

		    assert_eq!(encoder.bits2nuc::<u8>((encoder as u8).get_bits(6..8)), b'A');
		    assert_eq!(encoder.bits2nuc::<u8>((encoder as u8).get_bits(4..6)), b'C');
		    assert_eq!(encoder.bits2nuc::<u8>((encoder as u8).get_bits(2..4)), b'T');
		    assert_eq!(encoder.bits2nuc::<u8>((encoder as u8).get_bits(0..2)), b'G');
		)*
	    )
	}

        bits2one_base!(
            Naive::ACTG,
            Naive::ACGT,
            Naive::ATCG,
            Naive::ATGC,
            Naive::AGCT,
            Naive::AGTC,
            Naive::CATG,
            Naive::CAGT,
            Naive::CTAG,
            Naive::CTGA,
            Naive::CGAT,
            Naive::CGTA,
            Naive::TACG,
            Naive::TAGC,
            Naive::TCAG,
            Naive::TCGA,
            Naive::TGAC,
            Naive::TGCA,
            Naive::GACT,
            Naive::GATC,
            Naive::GCAT,
            Naive::GCTA,
            Naive::GTAC,
            Naive::GTCA
        );
    }

    #[test]
    fn comp_one_base_all_encoding() {
        macro_rules! comp_one_base {
	    ($($ty:expr), *) => (
		$(
		    let encoder = $ty;

		    assert_eq!(encoder.complement(encoder.nuc2bits::<u8>(b'A')), encoder.nuc2bits::<u8>(b'T'));
		    assert_eq!(encoder.complement(encoder.nuc2bits::<u8>(b'C')), encoder.nuc2bits::<u8>(b'G'));
		    assert_eq!(encoder.complement(encoder.nuc2bits::<u8>(b'T')), encoder.nuc2bits::<u8>(b'A'));
		    assert_eq!(encoder.complement(encoder.nuc2bits::<u8>(b'G')), encoder.nuc2bits::<u8>(b'C'));
		)*
	    )
	}

        comp_one_base!(
            Naive::ACTG,
            Naive::ACGT,
            Naive::ATCG,
            Naive::ATGC,
            Naive::AGCT,
            Naive::AGTC,
            Naive::CATG,
            Naive::CAGT,
            Naive::CTAG,
            Naive::CTGA,
            Naive::CGAT,
            Naive::CGTA,
            Naive::TACG,
            Naive::TAGC,
            Naive::TCAG,
            Naive::TCGA,
            Naive::TGAC,
            Naive::TGCA,
            Naive::GACT,
            Naive::GATC,
            Naive::GCAT,
            Naive::GCTA,
            Naive::GTAC,
            Naive::GTCA
        );
    }

    #[test]
    fn k15pu8() {
        let array = Naive::ACGT.encode(b"TAAGGATTCTAATCA");

        assert_eq!([131, 242, 13, 7], array);

        let table: Vec<u8> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);

        //One A more because encoder didn't know the size of kmer
        assert_eq!(Naive::ACGT.decode(array), b"TAAGGATTCTAATCAA");

        assert_eq!(
            Naive::ACGT.decode(Naive::ACGT.rev_comp(array)),
            b"TTGATTAGAATCCTTA"
        );
    }

    #[test]
    fn k15pu16() {
        let encoder = Naive::ACGT;

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([62083, 1805], array);

        let table: Vec<u16> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);

        //One more A because encoder didn't know the size of kmer
        assert_eq!(Naive::ACGT.decode(array), b"TAAGGATTCTAATCAA");

        assert_eq!(
            Naive::ACGT.decode(Naive::ACGT.rev_comp(array)),
            b"TTGATTAGAATCCTTA"
        );
    }

    #[test]
    fn k15pu32() {
        let encoder = Naive::ACGT;

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([118354563], array);

        let table: Vec<u32> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);

        //One more A because encoder didn't know the size of kmer
        assert_eq!(Naive::ACGT.decode(array), b"TAAGGATTCTAATCAA");

        assert_eq!(
            Naive::ACGT.decode(Naive::ACGT.rev_comp(array)),
            b"TTGATTAGAATCCTTA"
        );
    }

    #[test]
    fn k30pu32() {
        let encoder = Naive::ACGT;

        let array = encoder.encode(b"TAAGGATTCTAATCATAAGGATTCTAATCA");

        assert_eq!([3339580035, 29588640], array);

        let table: Vec<u32> = (0..30).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(
            table,
            vec![
                3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3,
                1, 0
            ]
        );

        //Two more A because encoder didn't know the size of kmer
        assert_eq!(
            Naive::ACGT.decode(array),
            b"TAAGGATTCTAATCATAAGGATTCTAATCAAA"
        );

        assert_eq!(
            Naive::ACGT.decode(Naive::ACGT.rev_comp(array)),
            b"TTTGATTAGAATCCTTATGATTAGAATCCTTA"
        );
    }

    #[test]
    fn k45pu64() {
        let encoder = Naive::ACGT;

        let array: [u64; kmer::word_for_k::<u64, 45>()] =
            encoder.encode(b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCA");

        assert_eq!([3585846758293238403, 7397160], array);

        let table: Vec<u64> = (0..45).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(
            table,
            vec![
                3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3,
                1, 0, 3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0
            ]
        );

        //19 more A because encoder didn't know the size of kmer
        assert_eq!(
            Naive::ACGT.decode(array),
            b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCAAAAAAAAAAAAAAAAAAAA"
        );

        assert_eq!(
            Naive::ACGT.decode(Naive::ACGT.rev_comp(array)),
            b"TTTTTTTTTTTTTTTTTTTTGATTAGAATCCTTATGATTAGAATCCTTATGATTAGAATCCTTA"
        );
    }

    #[test]
    fn k65pu128() {
        let encoder = Naive::ACGT;

        let array: [u128; kmer::word_for_k::<u128, 65>()] =
            encoder.encode(b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCAGGGGG");

        assert_eq!([226115275135941975929349834069397860995, 2], array);

        let table: Vec<u128> = (0..65).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(
            table,
            vec![
                3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3,
                1, 0, 3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0,
                0, 3, 1, 0, 2, 2, 2, 2, 2
            ]
        );

        // Many trailling A
        assert_eq!(Naive::ACGT.decode(array), b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCAGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        assert_eq!(Naive::ACGT.decode(Naive::ACGT.rev_comp(array)), b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCTGATTAGAATCCTTATGATTAGAATCCTTATGATTAGAATCCTTATGATTAGAATCCTTA");
    }
}
