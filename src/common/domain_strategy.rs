use franklin_crypto::bellman::{Engine, Field, PrimeField};

/// Padding prevents trivial collisions.
/// Each hash function nearly uses same padding strategies.
/// The only difference is that Rescue Prime requires no padding for
/// fixed length input. Rescue and Poseidon require same padding rule
/// for variable length input.
pub enum DomainStrategy<const RATE: usize> {
    // The capacity value is length x (^264 ) + (o − 1)
    // where o the output length. The padding consists of the field elements being 0.
    FixedLength,
    /// Padding is necessary for variable-length inputs, even if the input is already
    /// a multiple of the rate in length.
    // The capacity value is 2^64 + (o − 1) where o the output length.
    // The padding consists of one field element being 1,
    // and the remaining elements being 0
    VariableLength,
    // This is a variation of fixed length strategy. Only difference is value of capacity
    // element is being set to input length.    
    CustomFixedLength,
    CustomVariableLength,
    // No specialization and padding rule.
    NoPadding,
}

impl<const RATE: usize> DomainStrategy<RATE> {
    /// Computes capacity value for specialization and domain seperation.
    pub(crate) fn compute_capacity<E: Engine>(&self, input_len: usize) -> Option<E::Fr> {
        let mut repr = <E::Fr as PrimeField>::Repr::default();
        repr.as_mut()[1] = 1u64; // 2^64 corresponds second le limb
        let mut el = E::Fr::from_repr(repr).unwrap();

        let mut out_repr = <E::Fr as PrimeField>::Repr::default();
        out_repr.as_mut()[0] = (RATE - 1) as u64;
        let out_el = E::Fr::from_repr(repr).unwrap();

        match &self {
            Self::FixedLength => {
                // length * 2^64 + (o-1)
                // since we always use output length equals rate
                let length_as_fe = E::Fr::from_str(&input_len.to_string()).unwrap();
                el.mul_assign(&length_as_fe);
                el.add_assign(&out_el);

                Some(el)
            }
            Self::VariableLength => {
                // 2^64 + (o-1)
                el.add_assign(&out_el);

                Some(el)
            }
            Self::CustomFixedLength => {
                let mut repr = <E::Fr as PrimeField>::Repr::default();
                repr.as_mut()[0] = input_len as u64;

                E::Fr::from_repr(repr).ok()
            }
            Self::CustomVariableLength => None,
            _ => unimplemented!("unknown domain strategy"),
        }
    }
    /// Computes values for padding.
    pub(crate) fn generate_padding_values<E: Engine>(&self, input_len: usize) -> Vec<E::Fr> {
        assert!(input_len != 0, "empty input");
        let mut values_for_padding = vec![];
        match self {
            Self::FixedLength => {
                values_for_padding.resize(RATE - input_len, E::Fr::zero());

                values_for_padding
            }
            Self::VariableLength => {
                values_for_padding.push(E::Fr::one());
                while values_for_padding.len() % RATE != 0 {
                    values_for_padding.push(E::Fr::zero());
                }
                values_for_padding
            }

            Self::CustomFixedLength => {
                let mut cycle = input_len / RATE;

                if input_len % RATE != 0 {
                    cycle += 1;
                }

                let padding_len = cycle * RATE - input_len;

                for _ in 0..padding_len{
                    values_for_padding.push(E::Fr::zero());
                }
                
                values_for_padding
            }
            Self::CustomVariableLength => {
                assert!(input_len % RATE == 0, "input should be padded for variable length hashing");
                vec![]
            }
            _ => unimplemented!("unknown domain strategy"),
        }
    }
}
