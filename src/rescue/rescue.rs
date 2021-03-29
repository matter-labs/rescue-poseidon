use crate::common::{
    hash::generic_hash_with_padding, matrix::mmul_assign, padding::PaddingStrategy, sbox::sbox,
};
use crate::sponge::{SpongeMode, SpongeModes, SpongePermutation, SpongeState, StatefulSponge};
use crate::sponge_impl;
use crate::HasherParams;
use franklin_crypto::bellman::{Engine, Field};
/// The capacity value is length x (264 ) + (o − 1) where o the output length.
/// The padding consists of the field elements being 0.
pub fn rescue_hash_fixed_length<E: Engine, const S: usize, const R: usize>(
    input: &[E::Fr],
) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescueHasher<E, S, R>, S, R>(input, PaddingStrategy::FixedLength)
}

/// The capacity value is 264 + (o − 1) where o the output length.
/// The padding consists of one field element being 1, and the remaining elements being 0.
pub fn rescue_var_length<E: Engine, const S: usize, const R: usize>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescueHasher<E, S, R>, S, R>(
        input,
        PaddingStrategy::VariableLength,
    )
}

/// This is hasher with a custom strategy which basically sets value of capacity
pub fn rescue_hash<E: Engine, const S: usize, const R: usize>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, RescueHasher<E, S, R>, S, R>(input, PaddingStrategy::Custom)
}

#[derive(Debug, Clone)]
pub struct RescueHasher<E: Engine, const S: usize, const R: usize> {
    params: HasherParams<E, S, R>,
    state: [E::Fr; S],
    alpha: E::Fr,
    alpha_inv: E::Fr,
    sponge_mode: SpongeModes,
}

impl<E: Engine, const S: usize, const R: usize> Default for RescueHasher<E, S, R> {
    fn default() -> Self {
        let (params, alpha, alpha_inv) = super::params::rescue_params();
        Self {
            state: [E::Fr::zero(); S],
            params,
            alpha,
            alpha_inv: alpha_inv.expect("inverse of alpha"),
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}

impl<E: Engine, const S: usize, const R: usize> RescueHasher<E, S, R> {
    pub fn new_duplex() -> Self {
        let (params, alpha, alpha_inv) = super::params::rescue_params();
        Self {
            state: [E::Fr::zero(); S],
            alpha,
            alpha_inv: alpha_inv.expect("inverse of alpha"),
            params,
            sponge_mode: SpongeModes::Duplex(false),
        }
    }
}

// common parts of sponge
sponge_impl!(RescueHasher<E, S, R>);

impl<E: Engine, const S: usize, const R: usize> SpongePermutation<E> for RescueHasher<E, S, R> {
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
            mmul_assign::<E, S>(&self.params.mds_matrix, &mut self.state);

            // round constants
            self.state
                .iter_mut()
                .zip(self.params.constants_of_round(round + 1).iter())
                .for_each(|(s, c)| s.add_assign(c));
        }
    }
}
