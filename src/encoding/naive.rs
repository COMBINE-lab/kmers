//! Naive implementation of encoding
//! Not the best implementation but you can use any encoding

/* crate use */
use bit_field::BitArray as _;

/* Private function */

fn nuc2internal(nuc: u8) -> u8 {
    (nuc as u8 >> 1) & 0b11
}

const INTERNAL2NUC: [u8; 4] = [b'A', b'C', b'T', b'G'];

fn internal2nuc(internal: u8) -> u8 {
    INTERNAL2NUC[internal as usize]
}

const fn rev_encoding(encoding: u8) -> u8 {
    let mut rev = 0;

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

    pub(crate) fn bits2nuc<P>(&self, bits: P) -> u8
    where
        P: crate::utils::Data,
    {
        // I forget how this work but it's works
        let reverse = rev_encoding(*self as u8);
        internal2nuc((reverse >> (6 - (bits.to_u8() & 0b11) * 2)) & 0b11)
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

    /// Convert an array of two bits data in
    fn decode(&self, array: [P; B]) -> Vec<u8> {
        let mut seq = Vec::with_capacity(B * 4);

        for idx in 0..P::BIT_LENGTH * B / 2 {
            let value = array.get_bits(idx * 2..=idx * 2 + 1);
            seq.push(self.bits2nuc(value));
        }

        seq
    }

    fn rev_comp(&self, _array: [P; B]) -> [P; B] {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::encoding::Encoding as _;

    use crate::kmer;

    #[test]
    fn one_base_all_encoding() {
        let encoder = Naive::ACTG;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b11);

        let encoder = Naive::ACGT;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b10);

        let encoder = Naive::ATCG;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b11);

        let encoder = Naive::ATGC;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b10);

        let encoder = Naive::AGCT;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b01);

        let encoder = Naive::AGTC;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b01);

        let encoder = Naive::CATG;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b11);

        let encoder = Naive::CAGT;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b10);

        let encoder = Naive::CTAG;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b11);

        let encoder = Naive::CTGA;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b10);

        let encoder = Naive::CGAT;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b01);

        let encoder = Naive::CGTA;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b01);

        let encoder = Naive::TACG;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b11);

        let encoder = Naive::TAGC;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b10);

        let encoder = Naive::TCAG;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b11);

        let encoder = Naive::TCGA;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b10);

        let encoder = Naive::TGAC;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b01);

        let encoder = Naive::TGCA;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b01);

        let encoder = Naive::GACT;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b00);

        let encoder = Naive::GATC;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b00);

        let encoder = Naive::GCAT;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b00);

        let encoder = Naive::GCTA;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b00);

        let encoder = Naive::GTAC;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b00);

        let encoder = Naive::GTCA;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b11);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b00);
    }

    #[test]
    fn one_base_all_decoding() {
        let encoder = Naive::ACTG;

        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'G');

        let encoder = Naive::ACGT;

        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'G');

        let encoder = Naive::ATCG;

        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'G');

        let encoder = Naive::ATGC;

        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'G');

        let encoder = Naive::AGCT;

        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'G');

        let encoder = Naive::AGTC;

        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'G');

        let encoder = Naive::CATG;

        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'G');

        let encoder = Naive::CAGT;

        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'G');

        let encoder = Naive::CTAG;

        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'G');

        let encoder = Naive::CTGA;

        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'G');

        let encoder = Naive::CGAT;

        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'G');

        let encoder = Naive::CGTA;

        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'G');

        let encoder = Naive::TACG;

        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'G');

        let encoder = Naive::TAGC;

        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'G');

        let encoder = Naive::TCAG;

        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'G');

        let encoder = Naive::TCGA;

        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'G');

        let encoder = Naive::TGAC;

        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'G');

        let encoder = Naive::TGCA;

        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'G');

        let encoder = Naive::GACT;

        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'G');

        let encoder = Naive::GATC;

        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'G');

        let encoder = Naive::GCAT;

        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'G');

        let encoder = Naive::GCTA;

        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'G');

        let encoder = Naive::GTAC;

        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'G');

        let encoder = Naive::GTCA;

        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'G');
    }

    #[test]
    fn k15pu8() {
        let array = Naive::ACGT.encode(b"TAAGGATTCTAATCA");

        assert_eq!([131, 242, 13, 7], array);

        let table: Vec<u8> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);

        assert_eq!(Naive::ACGT.decode(array), b"TAAGGATTCTAATCAA"); //One A more because encoder didn't know the size of kmer
    }

    #[test]
    fn k15pu16() {
        let encoder = Naive::ACGT;

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([62083, 1805], array);

        let table: Vec<u16> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);
    }

    #[test]
    fn k15pu32() {
        let encoder = Naive::ACGT;

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([118354563], array);

        let table: Vec<u32> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![3, 0, 0, 2, 2, 0, 3, 3, 1, 3, 0, 0, 3, 1, 0]);
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
    }
}
