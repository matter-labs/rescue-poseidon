use crate::{common::matrix::mmul_assign};
use crate::common::{hash::generic_hash_with_padding, padding::PaddingStrategy, sbox::sbox};
use crate::sponge::{SpongePermutation, SpongeState, StatefulSponge, SpongeMode, SpongeModes};
use crate::sponge_impl;
use crate::HasherParams;
use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;

/// The capacity value is length x (264 ) + (o − 1) where o the output length.
/// The padding consists of the field elements being 0.
pub fn rescue_prime_fixed_length<E: Engine, const S: usize, const R: usize>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescuePrimeHasher<E, S, R>, S, R>(
        input,
        PaddingStrategy::FixedLength,
    )
}

/// The capacity value is 264 + (o − 1) where o the output length.
/// The padding consists of one field element being 1, and the remaining elements being 0.
pub fn rescue_prime_var_length<E: Engine, const S: usize, const R: usize>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescuePrimeHasher<E, S, R>,S, R>(
        input,
        PaddingStrategy::VariableLength,
    )
}

/// This is hasher with a custom strategy which basically sets value of capacity
pub fn rescue_prime_hash<E: Engine, const S: usize, const R: usize>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescuePrimeHasher<E, S, R>,S, R>(
        input,
        PaddingStrategy::Custom,
    )
}

#[derive(Debug, Clone)]
pub struct RescuePrimeHasher<E: Engine, const S: usize, const R: usize> {
    params: HasherParams<E, R, S>,
    state: [E::Fr; S],
    alpha: E::Fr,
    alpha_inv: E::Fr,
    sponge_mode: SpongeModes,
}

impl<E: Engine, const S: usize, const R: usize> Default for RescuePrimeHasher<E, S, R> {
    fn default() -> Self {
        let (params, alpha, alpha_inv) = crate::rescue_prime::params::rescue_prime_params();
        Self {
            state: [E::Fr::zero(); S],
            params,
            alpha,
            alpha_inv,
            // TODO
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}

// common parts of sponge
sponge_impl!(RescuePrimeHasher<E, S, R>);
// sponge_impl!(RescuePrimeHasher<E>, super::params::rescue_prime_params);

impl<E: Engine, const S: usize, const R: usize> SpongePermutation<E> for RescuePrimeHasher<E,S, R> {
    fn permutation(&mut self) {
        for round in 0..self.params.full_rounds {
            // sbox alpha
            sbox::<E>(self.alpha, &mut self.state);
            // mds
            mmul_assign::<E, S>(&self.params.mds_matrix, &mut self.state);

            // round constants
            self.state
                .iter_mut()
                .zip(self.params.constants_of_round(round).iter())
                .for_each(|(s, c)| s.add_assign(c));

            // sbox alpha inv
            sbox::<E>(self.alpha_inv, &mut self.state);

            // mds
            mmul_assign::<E, S>(&self.params.mds_matrix, &mut self.state);

            // round constants
            self.state
                .iter_mut()
                .zip(self.params.constants_of_round(round + 1).iter())
                .for_each(|(s, c)| s.add_assign(c));
        }
    }
}
