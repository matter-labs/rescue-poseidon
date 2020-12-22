use crate::{common::matrix::mmul_assign};
use crate::common::{hash::generic_hash_with_padding, padding::PaddingStrategy, sbox::sbox};
use crate::sponge::{SpongeParams, SpongePermutation, SpongeState, StatefulSponge, SpongeMode, SpongeModes};
use crate::sponge_impl;
use crate::HasherParams;
use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;

/// The capacity value is length x (264 ) + (o − 1) where o the output length.
/// The padding consists of the field elements being 0.
pub fn rescue_prime_fixed_length<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescuePrimeHasher<E>>(
        input,
        PaddingStrategy::FixedLength(input.len()),
    )
}

/// The capacity value is 264 + (o − 1) where o the output length.
/// The padding consists of one field element being 1, and the remaining elements being 0.
pub fn rescue_prime_var_length<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescuePrimeHasher<E>>(
        input,
        PaddingStrategy::VariableLength(input.len()),
    )
}

/// This is hasher with a custom strategy which basically sets value of capacity
pub fn rescue_prime_hash<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescuePrimeHasher<E>>(
        input,
        PaddingStrategy::Custom(input.len()),
    )
}

#[derive(Debug, Clone)]
pub struct RescuePrimeHasher<E: Engine> {
    params: HasherParams<E>,
    state: Vec<E::Fr>,
    tmp_storage: Vec<E::Fr>,
    alpha: E::Fr,
    alpha_inv: E::Fr,
    sponge_mode: SpongeModes<E>,
}

impl<E: Engine> Default for RescuePrimeHasher<E> {
    fn default() -> Self {
        let (params, alpha, alpha_inv) = crate::rescue_prime::params::rescue_prime_params();
        Self {
            state: vec![E::Fr::zero(); params.state_width],
            tmp_storage: Vec::with_capacity(params.rate),
            params,
            alpha,
            alpha_inv,
            sponge_mode: SpongeModes::Standard,
        }
    }
}

// common parts of sponge
sponge_impl!(RescuePrimeHasher<E>);
// sponge_impl!(RescuePrimeHasher<E>, super::params::rescue_prime_params);

impl<E: Engine> SpongePermutation<E> for RescuePrimeHasher<E> {
    fn permutation(&mut self) {
        for round in 0..self.params.full_rounds {
            // sbox alpha
            sbox::<E>(self.alpha, &mut self.state);
            // mds
            mmul_assign::<E>(&self.params.mds_matrix, &mut self.state);

            // round constants
            self.state
                .iter_mut()
                .zip(self.params.constants_of_round(round).iter())
                .for_each(|(s, c)| s.add_assign(c));

            // sbox alpha inv
            sbox::<E>(self.alpha_inv, &mut self.state);

            // mds
            mmul_assign::<E>(&self.params.mds_matrix, &mut self.state);

            // round constants
            self.state
                .iter_mut()
                .zip(self.params.constants_of_round(round + 1).iter())
                .for_each(|(s, c)| s.add_assign(c));
        }
    }
}
