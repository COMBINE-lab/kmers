//! XOR10 implementation of encoding
//! A fast encoder but A -> 00, C -> 01, T -> 10, G -> 11

/* crate use */
use bit_field::BitArray as _;

/// Lookup table usefull to convert internal encoding in ASCII
const BITS2NUC: [u8; 4] = [b'A', b'C', b'T', b'G'];

pub struct Xor10;

impl Xor10 {
    /// Convert nucleotide in encoding corresponding 2 bits
    #[inline]
    pub(crate) fn nuc2bits<P>(&self, nuc: u8) -> P
    where
        P: crate::utils::Data,
    {
        P::from((nuc >> 1) & 0b11)
    }

    /// Convert nucleotide encode on 2 bits in ASCII
    #[inline]
    pub(crate) fn bits2nuc<P>(&self, bits: P) -> u8
    where
        P: crate::utils::Data,
    {
        BITS2NUC[bits.to_u8() as usize]
    }

    /// Get the complement of a nucleotide encode in 2 bits
    #[inline]
    pub(crate) fn complement<P>(&self, bits: P) -> P
    where
        P: crate::utils::Data,
    {
        P::from(bits.to_u8() ^ 0b10)
    }
}

impl<P, const B: usize> super::Encoding<P, B> for Xor10
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

    fn rev_comp<const K: usize>(&self, mut array: [P; B]) -> [P; B] {
        // This could probably be improve natir/cocktail have a nicer implementation for u64
        let mut i = 0;
        let mut j = K * 2 - 2;

        while i <= j {
            let comp_i = self.complement(array.get_bits(i..i + 2));
            let comp_j = self.complement(array.get_bits(j..j + 2));

            array.set_bits(i..i + 2, comp_j);
            array.set_bits(j..j + 2, comp_i);

            i += 2;
            j -= 2;
        }

        array
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /* Project use */
    use crate::kmer;

    use crate::encoding::Encoding as _;

    #[test]
    fn one_base_encoding() {
        let encoder = Xor10;

        assert_eq!(encoder.nuc2bits::<u8>(b'A'), 0b00);
        assert_eq!(encoder.nuc2bits::<u8>(b'C'), 0b01);
        assert_eq!(encoder.nuc2bits::<u8>(b'T'), 0b10);
        assert_eq!(encoder.nuc2bits::<u8>(b'G'), 0b11);
    }

    #[test]
    fn one_base_decoding() {
        let encoder = Xor10;

        assert_eq!(encoder.bits2nuc::<u8>(0b00), b'A');
        assert_eq!(encoder.bits2nuc::<u8>(0b01), b'C');
        assert_eq!(encoder.bits2nuc::<u8>(0b10), b'T');
        assert_eq!(encoder.bits2nuc::<u8>(0b11), b'G');
    }

    #[test]
    fn one_base_comp() {
        let encoder = Xor10;

        assert_eq!(
            encoder.complement(encoder.nuc2bits::<u8>(b'A')),
            encoder.nuc2bits::<u8>(b'T')
        );
        assert_eq!(
            encoder.complement(encoder.nuc2bits::<u8>(b'C')),
            encoder.nuc2bits::<u8>(b'G')
        );
        assert_eq!(
            encoder.complement(encoder.nuc2bits::<u8>(b'T')),
            encoder.nuc2bits::<u8>(b'A')
        );
        assert_eq!(
            encoder.complement(encoder.nuc2bits::<u8>(b'G')),
            encoder.nuc2bits::<u8>(b'C')
        );
    }

    #[test]
    fn k15pu8() {
        let array = Xor10.encode(b"TAAGGATTCTAATCA");

        assert_eq!([194, 163, 9, 6], array);

        let table: Vec<u8> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0]);

        //One A more because encoder didn't know the size of kmer
        assert_eq!(Xor10.decode(array), b"TAAGGATTCTAATCAA");

        assert_eq!(
            Xor10.decode(Xor10.rev_comp::<15>(array)),
            b"TGATTAGAATCCTTAA"
        );
    }

    #[test]
    fn k15pu16() {
        let encoder = Xor10;

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([41922, 1545], array);

        let table: Vec<u16> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0]);

        //One more A because encoder didn't know the size of kmer
        assert_eq!(Xor10.decode(array), b"TAAGGATTCTAATCAA");

        assert_eq!(
            Xor10.decode(Xor10.rev_comp::<15>(array)),
            b"TGATTAGAATCCTTAA"
        );
    }

    #[test]
    fn k15pu32() {
        let encoder = Xor10;

        let array = encoder.encode(b"TAAGGATTCTAATCA");

        assert_eq!([101295042], array);

        let table: Vec<u32> = (0..15).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(table, vec![2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0]);

        //One more A because encoder didn't know the size of kmer
        assert_eq!(Xor10.decode(array), b"TAAGGATTCTAATCAA");

        assert_eq!(
            Xor10.decode(Xor10.rev_comp::<15>(array)),
            b"TGATTAGAATCCTTAA"
        );
    }

    #[test]
    fn k30pu32() {
        let encoder = Xor10;

        let array = encoder.encode(b"TAAGGATTCTAATCATAAGGATTCTAATCA");

        assert_eq!([2248778690, 25323760], array);

        let table: Vec<u32> = (0..30).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(
            table,
            vec![
                2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0, 2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2,
                1, 0
            ]
        );

        //Two more A because encoder didn't know the size of kmer
        assert_eq!(Xor10.decode(array), b"TAAGGATTCTAATCATAAGGATTCTAATCAAA");

        assert_eq!(
            Xor10.decode(Xor10.rev_comp::<30>(array)),
            b"TGATTAGAATCCTTATGATTAGAATCCTTAAA"
        );
    }

    #[test]
    fn k45pu64() {
        let encoder = Xor10;

        let array: [u64; kmer::word_for_k::<u64, 45>()] =
            encoder.encode(b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCA");

        assert_eq!([2414607732474225602, 6330940], array);

        let table: Vec<u64> = (0..45).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(
            table,
            vec![
                2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0, 2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2,
                1, 0, 2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0
            ]
        );

        //19 more A because encoder didn't know the size of kmer
        assert_eq!(
            Xor10.decode(array),
            b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCAAAAAAAAAAAAAAAAAAAA"
        );

        assert_eq!(
            Xor10.decode(Xor10.rev_comp::<45>(array)),
            b"TGATTAGAATCCTTATGATTAGAATCCTTATGATTAGAATCCTTAAAAAAAAAAAAAAAAAAAA"
        );
    }

    #[test]
    fn k65pu128() {
        let encoder = Xor10;

        let array: [u128; kmer::word_for_k::<u128, 65>()] =
            encoder.encode(b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCAGGGGG");

        assert_eq!([339078536113543227067743297186703188930, 3], array);

        let table: Vec<u128> = (0..65).map(|i| array.get_bits(i * 2..=i * 2 + 1)).collect();

        assert_eq!(
            table,
            vec![
                2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0, 2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2,
                1, 0, 2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0, 0, 2, 1, 0, 2, 0, 0, 3, 3, 0, 2, 2, 1, 2, 0,
                0, 2, 1, 0, 3, 3, 3, 3, 3
            ]
        );

        // Many trailling A
        assert_eq!(Xor10.decode(array), b"TAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCATAAGGATTCTAATCAGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        assert_eq!(Xor10.decode(Xor10.rev_comp::<65>(array)), b"CCCCCTGATTAGAATCCTTATGATTAGAATCCTTATGATTAGAATCCTTATGATTAGAATCCTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    }
}
