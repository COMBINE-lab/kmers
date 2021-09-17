/* crate use */
use bit_field::BitArray;

/* project use */
use super::*;

/// Enumeration of all possible encoding:
/// - the first nucleotide is equal to 00
/// - the second nucleotide is equal to 01
/// - the third nucleotide is equal to 10
/// - the last nucleotide is equal to 11
#[derive(Clone, Copy, Debug)]
pub enum Encoding {
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

const fn nuc2internal(nuc: u8) -> u8 {
    (nuc as u8 >> 1) & 0b11
}

pub struct Naive<P> {
    encoding: Encoding,
    _marker: std::marker::PhantomData<P>,
}

impl<P> Naive<P>
where
    P: std::convert::From<u8>,
{
    pub fn new(encoding: Encoding) -> Self {
        Self {
            encoding,
            _marker: std::marker::PhantomData,
        }
    }

    /// Convert nucleotide in encoding corresponding 2 bits  
    pub(crate) fn nuc2bits(&self, nuc: u8) -> P {
        let index = 6 - nuc2internal(nuc) * 2;
	
        P::from((self.encoding as u8 >> index) & 0b11)
    }
}

impl<P, const B: usize> Encoder<P, B> for Naive<P>
where
    P: std::convert::From<u8> + bit_field::BitField + std::marker::Copy,
{
    fn encode(&self, seq: &[u8]) -> [P; B] {
        let mut array: [P; B] = unsafe { [std::mem::zeroed(); B] };

        for (idx, nuc) in seq.iter().enumerate() {
            array.set_bits(idx * 2..=idx * 2 + 1, self.nuc2bits(*nuc));
        }

        array
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bit_field::BitArray;

    use crate::kmer;

    #[test]
    fn one_base_encoding() {
        let encoder = Naive::<u8>::new(Encoding::ACTG);

        assert_eq!(encoder.nuc2bits(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits(b'G'), 0b11);

        let encoder = Naive::<u8>::new(Encoding::ACGT);

        assert_eq!(encoder.nuc2bits(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits(b'G'), 0b10);

        let encoder = Naive::<u8>::new(Encoding::ATCG);

        assert_eq!(encoder.nuc2bits(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits(b'G'), 0b11);

        let encoder = Naive::<u8>::new(Encoding::ATGC);

        assert_eq!(encoder.nuc2bits(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits(b'G'), 0b10);

        let encoder = Naive::<u8>::new(Encoding::AGCT);

        assert_eq!(encoder.nuc2bits(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits(b'G'), 0b01);

        let encoder = Naive::<u8>::new(Encoding::AGTC);

        assert_eq!(encoder.nuc2bits(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits(b'G'), 0b01);

        let encoder = Naive::<u8>::new(Encoding::CATG);

        assert_eq!(encoder.nuc2bits(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits(b'G'), 0b11);

        let encoder = Naive::<u8>::new(Encoding::CAGT);

        assert_eq!(encoder.nuc2bits(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits(b'G'), 0b10);

        let encoder = Naive::<u8>::new(Encoding::CTAG);

        assert_eq!(encoder.nuc2bits(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits(b'G'), 0b11);

        let encoder = Naive::<u8>::new(Encoding::CTGA);

        assert_eq!(encoder.nuc2bits(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits(b'G'), 0b10);

        let encoder = Naive::<u8>::new(Encoding::CGAT);

        assert_eq!(encoder.nuc2bits(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits(b'G'), 0b01);

        let encoder = Naive::<u8>::new(Encoding::CGTA);

        assert_eq!(encoder.nuc2bits(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits(b'G'), 0b01);

        let encoder = Naive::<u8>::new(Encoding::TACG);

        assert_eq!(encoder.nuc2bits(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits(b'G'), 0b11);

        let encoder = Naive::<u8>::new(Encoding::TAGC);

        assert_eq!(encoder.nuc2bits(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits(b'G'), 0b10);

        let encoder = Naive::<u8>::new(Encoding::TCAG);

        assert_eq!(encoder.nuc2bits(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits(b'G'), 0b11);

        let encoder = Naive::<u8>::new(Encoding::TCGA);

        assert_eq!(encoder.nuc2bits(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits(b'G'), 0b10);

        let encoder = Naive::<u8>::new(Encoding::TGAC);

        assert_eq!(encoder.nuc2bits(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits(b'G'), 0b01);

        let encoder = Naive::<u8>::new(Encoding::TGCA);

        assert_eq!(encoder.nuc2bits(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits(b'G'), 0b01);

        let encoder = Naive::<u8>::new(Encoding::GACT);

        assert_eq!(encoder.nuc2bits(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits(b'G'), 0b00);

        let encoder = Naive::<u8>::new(Encoding::GATC);

        assert_eq!(encoder.nuc2bits(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits(b'G'), 0b00);

        let encoder = Naive::<u8>::new(Encoding::GCAT);

        assert_eq!(encoder.nuc2bits(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits(b'G'), 0b00);

        let encoder = Naive::<u8>::new(Encoding::GCTA);

        assert_eq!(encoder.nuc2bits(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits(b'G'), 0b00);

        let encoder = Naive::<u8>::new(Encoding::GTAC);

        assert_eq!(encoder.nuc2bits(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits(b'G'), 0b00);

        let encoder = Naive::<u8>::new(Encoding::GTCA);

        assert_eq!(encoder.nuc2bits(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits(b'G'), 0b00);
    }

    #[test]
    fn k15pu8() {
        let encoder = Naive::<u8>::new(Encoding::ACGT);

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([131, 242, 13, 7], array);

        let table: Vec<u8> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);
    }

    #[test]
    fn k15pu16() {
        let encoder = Naive::<u16>::new(Encoding::ACGT);

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([62083, 1805], array);

        let table: Vec<u16> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);
    }

    #[test]
    fn k15pu32() {
        let encoder = Naive::<u32>::new(Encoding::ACGT);

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([118354563], array);

        let table: Vec<u32> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);
    }

    #[test]
    fn k30pu32() {
        let encoder = Naive::<u32>::new(Encoding::ACGT);

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
    }

    #[test]
    fn k45pu64() {
        let encoder = Naive::<u64>::new(Encoding::ACGT);

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
    }

    #[test]
    fn k65pu128() {
        let encoder = Naive::<u128>::new(Encoding::ACGT);

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
    }
}
