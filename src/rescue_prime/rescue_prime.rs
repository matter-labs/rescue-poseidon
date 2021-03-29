use crate::{common::matrix::mmul_assign};
use crate::common::{hash::generic_hash_with_padding, domain_strategy::DomainStrategy, sbox::sbox};
use crate::sponge::{SpongePermutation, SpongeState, StatefulSponge, SpongeMode, SpongeModes};
use crate::sponge_impl;
use crate::common::params::HasherParams;
use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;
use std::convert::TryInto;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter. Uses state-width=3 and rate=2.
pub fn rescue_prime_hash<E: Engine, const L: usize>(input: &[E::Fr; L]) -> [E::Fr; 2] {
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;
    
    rescue_prime_generic_fixed_length::<E, STATE_WIDTH, RATE, L>(input)
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule. Uses state-width=3 and rate=2.
pub fn rescue_prime_hash_var_length<E: Engine>(input: &[E::Fr]) -> [E::Fr; 2] {
    // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;

    rescue_prime_generic_var_length::<E, STATE_WIDTH, RATE>(input)
}

pub(crate) fn rescue_prime_generic_fixed_length<
    E: Engine,
    const STATE_WIDTH: usize,
    const RATE: usize,
    const LENGTH: usize,
>(
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    let result =
        generic_hash_with_padding::<E, RescuePrimeHasher<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            input,
            DomainStrategy::CustomFixedLength,
        );

    result.try_into().expect("fixed length array")
}

pub(crate) fn rescue_prime_generic_var_length<E: Engine, const STATE_WIDTH: usize, const RATE: usize>(
    input: &[E::Fr],
) -> [E::Fr; RATE] {
    let result =
        generic_hash_with_padding::<E, RescuePrimeHasher<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            input,
            DomainStrategy::CustomVariableLength,
        );

    result.try_into().expect("fixed length array")
}

#[derive(Debug, Clone)]
pub struct RescuePrimeHasher<E: Engine, const S: usize, const R: usize> {
    params: HasherParams<E, S, R>,
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
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}
impl<E: Engine, const S: usize, const R: usize> RescuePrimeHasher<E, S, R> {
    fn new_duplex() -> Self {
        let (params, alpha, alpha_inv) = crate::rescue_prime::params::rescue_prime_params();
        Self {
            state: [E::Fr::zero(); S],
            params,
            alpha,
            alpha_inv,
            sponge_mode: SpongeModes::Duplex(false),
        }
    }
}

// common parts of sponge
sponge_impl!(RescuePrimeHasher<E, S, R>);

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
