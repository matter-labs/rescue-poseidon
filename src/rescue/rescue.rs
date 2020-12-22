use crate::{common::{
    hash::generic_hash_with_padding, matrix::mmul_assign, padding::PaddingStrategy, sbox::sbox,
}};
use crate::sponge::{SpongeParams, SpongePermutation, SpongeState, SpongeMode, SpongeModes, StatefulSponge};
use crate::sponge_impl;
use crate::HasherParams;
use franklin_crypto::bellman::{Engine, Field};
/// The capacity value is length x (264 ) + (o − 1) where o the output length.
/// The padding consists of the field elements being 0.
pub fn rescue_fixed_length<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescueHasher<E>>(
        input,
        PaddingStrategy::FixedLength(input.len()),
    )
}

/// The capacity value is 264 + (o − 1) where o the output length.
/// The padding consists of one field element being 1, and the remaining elements being 0.
pub fn rescue_var_length<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescueHasher<E>>(
        input,
        PaddingStrategy::VariableLength(input.len()),
    )
}

/// This is hasher with a custom strategy which basically sets value of capacity
pub fn rescue_hash<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescueHasher<E>>(input, PaddingStrategy::Custom(input.len()))
}

#[derive(Debug, Clone)]
pub struct RescueHasher<E: Engine> {
    params: HasherParams<E>,
    state: Vec<E::Fr>,
    tmp_storage: Vec<E::Fr>,
    alpha: E::Fr,
    alpha_inv: E::Fr,
    sponge_mode: SpongeModes<E>,
}

impl<E: Engine> Default for RescueHasher<E> {
    fn default() -> Self {
        let (params, alpha, alpha_inv) = super::params::rescue_params();
        Self {
            state: vec![E::Fr::zero(); params.state_width],
            tmp_storage: Vec::with_capacity(params.rate),
            params,
            alpha,
            alpha_inv: alpha_inv.expect("inverse of alpha"),
            sponge_mode: SpongeModes::Standard,
        }
    }
}

impl<E: Engine> RescueHasher<E> {
    pub fn new(sponge_mode: SpongeModes<E>) -> Self {
        let (params, alpha, alpha_inv) = super::params::rescue_params();
        Self {
            state: vec![E::Fr::zero(); params.state_width],
            tmp_storage: Vec::with_capacity(params.rate),
            alpha,
            alpha_inv: alpha_inv.expect("inverse of alpha"),
            params,
            sponge_mode,
        }
    }
}

// common parts of sponge
sponge_impl!(RescueHasher<E>);
// sponge_impl!(RescueHasher<E>, super::params::rescue_params);

impl<E: Engine> SpongePermutation<E> for RescueHasher<E> {
    fn permutation(&mut self) {
        // round constants for first step
        self.state
            .iter_mut()
            .zip(self.params.constants_of_round(0).iter())
            .for_each(|(s, c)| s.add_assign(c));

        for round in 0..2 * self.params.full_rounds {
            // sbox
            if round & 1 == 0 {
                sbox::<E>(self.alpha_inv, &mut self.state);
            } else {
                sbox::<E>(self.alpha, &mut self.state);
            }

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
