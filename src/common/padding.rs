use franklin_crypto::bellman::{Engine, Field, PrimeField};

/// Padding prevents trivial collisions.
/// Each hash function nearly uses same padding strategies.
/// The only difference is that Rescue Prime requires no padding for
/// fixed length input. Rescue and Poseidon require same padding rule
/// for variable length input.
pub enum PaddingStrategy {
    // The capacity value is length x (^264 ) + (o − 1)
    // where o the output length. The padding consists of the field elements being 0.
    FixedLength,
    /// Padding is necessary for variable-length inputs, even if the input is already
    /// a multiple of the rate in length.
    // The capacity value is 2^64 + (o − 1) where o the output length.
    // The padding consists of one field element being 1,
    // and the remaining elements being 0
    VariableLength,
    // zksync uses a custom specialization which basically sets value of capacity element
    // to te input length. The only difference from variable length strategy, this is applied 
    // when input length is not multiple of rate param.
    Custom,
    // No specialization and padding rule.
    NoPadding,
}

impl PaddingStrategy {
    /// Computes capacity value for specialization and domain seperation.
    pub(crate) fn compute_capacity<E: Engine>(&self, input_len: usize, rate: usize) -> Option<E::Fr> {
        let mut repr = <E::Fr as PrimeField>::Repr::default();
        repr.as_mut()[1] = 1u64; // 2^64 corresponds second le limb
        let mut el = E::Fr::from_repr(repr).unwrap();

        let mut out_repr = <E::Fr as PrimeField>::Repr::default();
        out_repr.as_mut()[0] = (rate - 1) as u64;
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
            Self::Custom => {
                let mut repr = <E::Fr as PrimeField>::Repr::default();
                repr.as_mut()[0] = input_len as u64;

                E::Fr::from_repr(repr).ok()
            }
            _ => None,
        }
    }
    /// Computes values for padding.
    pub(crate) fn generate_padding_values<E: Engine>(&self, input_len: usize, rate: usize) -> Vec<E::Fr> {
        let mut values_for_padding = vec![];
        match &self {
            Self::FixedLength => {
                values_for_padding.resize(rate - input_len, E::Fr::zero());

                values_for_padding
            }
            Self::VariableLength => {
                values_for_padding.push(E::Fr::one());
                while values_for_padding.len() % rate != 0 {
                    values_for_padding.push(E::Fr::zero());
                }
                values_for_padding
            }
            Self::Custom => {
                if rate - input_len > 0 {
                    values_for_padding.push(E::Fr::one());
                }
                while values_for_padding.len() % rate != 0 {
                    values_for_padding.push(E::Fr::zero());
                }

                values_for_padding
            }
            Self::NoPadding => values_for_padding,
        }
    }
}
