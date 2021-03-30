use super::domain_strategy::DomainStrategy;
use crate::sponge::StatefulSponge;
use franklin_crypto::bellman::Engine;

// This function specializes sponge construction with respect to domain strategy.
pub(crate) fn generic_hash_with_padding<
    E: Engine,
    Sponge: StatefulSponge<E, S, R>,
    const S: usize,
    const R: usize,
>(
    input: &[E::Fr],
    domain_strategy: DomainStrategy<R>,
) -> Vec<E::Fr> {
    let mut sponge = Sponge::default();

    let capacity_value = domain_strategy.compute_capacity::<E>(input.len());

    let padding_values = domain_strategy.generate_padding_values::<E>(input.len());

    let input_with_padding = input
        .iter()
        .chain(&padding_values[..])
        .cloned()
        .collect::<Vec<E::Fr>>();

    sponge.specialize(capacity_value);
    sponge.absorb(&input_with_padding);
    sponge.squeeze(None)
}

#[macro_export]
macro_rules! hash_imp {
    ($hasher_name:ty) => {
        /// Receives inputs whose length `known` prior(fixed-length).
        /// Also uses custom domain strategy which basically sets value of capacity element to
        /// length of input and applies a padding rule which makes input size equals to multiple of
        /// rate parameter. Uses state-width=3 and rate=2.
        pub fn $hasher_name<E: Engine, const L: usize>(input: &[E::Fr; L]) -> [E::Fr; 2] {
            const STATE_WIDTH: usize = 3;
            const RATE: usize = 2;

            poseidon_generic_fixed_length::<E, STATE_WIDTH, RATE, L>(input)
        }

        /// Receives inputs whose length `unknown` prior (variable-length).
        /// Also uses custom domain strategy which does not touch to value of capacity element
        /// and does not apply any padding rule. Uses state-width=3 and rate=2.
        pub fn poseidon_hash_var_length<E: Engine>(input: &[E::Fr]) -> [E::Fr; 2] {
            // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
            const STATE_WIDTH: usize = 3;
            const RATE: usize = 2;

            poseidon_generic_var_length::<E, STATE_WIDTH, RATE>(input)
        }
    };
}
