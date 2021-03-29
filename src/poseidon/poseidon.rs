use crate::common::domain_strategy::DomainStrategy;
use crate::common::matrix::mmul_assign;
use crate::sponge_impl;
use crate::{common::hash::generic_hash_with_padding, common::params::HasherParams};
use crate::{
    common::sbox::sbox,
    sponge::{SpongeMode, SpongeModes, SpongePermutation, SpongeState, StatefulSponge},
};
use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;
use std::convert::TryInto;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter. Uses state-width=3 and rate=2.
pub fn poseidon_hash<E: Engine, const L: usize>(input: &[E::Fr; L]) -> [E::Fr; 2] {
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

pub(crate) fn poseidon_generic_fixed_length<
    E: Engine,
    const STATE_WIDTH: usize,
    const RATE: usize,
    const LENGTH: usize,
>(
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    let result =
        generic_hash_with_padding::<E, PoseidonHasher<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            input,
            DomainStrategy::CustomFixedLength,
        );

    result.try_into().expect("fixed length array")
}

pub(crate) fn poseidon_generic_var_length<
    E: Engine,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    input: &[E::Fr],
) -> [E::Fr; RATE] {
    let result =
        generic_hash_with_padding::<E, PoseidonHasher<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            input,
            DomainStrategy::CustomVariableLength,
        );

    result.try_into().expect("fixed length array")
}

#[derive(Clone, Debug)]
pub struct PoseidonHasher<E: Engine, const S: usize, const R: usize> {
    params: HasherParams<E, S, R>,
    state: [E::Fr; S],
    alpha: E::Fr,
    optimized_round_constants: Vec<[E::Fr; S]>,
    optimized_mds_matrixes: ([[E::Fr; S]; S], Vec<[[E::Fr; S]; S]>),
    sponge_mode: SpongeModes,
}

impl<E: Engine, const S: usize, const R: usize> Default for PoseidonHasher<E, S, R> {
    fn default() -> Self {
        let (params, alpha, optimized_round_constants, optimized_mds_matrixes) =
            super::params::poseidon_light_params();
        Self {
            state: [E::Fr::zero(); S],
            alpha,
            params,
            optimized_round_constants,
            optimized_mds_matrixes,
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}
impl<E: Engine, const S: usize, const R: usize> PoseidonHasher<E, S, R> {
    pub fn new_duplex() -> Self {
        let (params, alpha, optimized_round_constants, optimized_mds_matrixes) =
            super::params::poseidon_light_params();
        Self {
            state: [E::Fr::zero(); S],
            alpha,
            params,
            optimized_round_constants,
            optimized_mds_matrixes,
            sponge_mode: SpongeModes::Duplex(false),
        }
    }
}

// Implementation of common parts
sponge_impl!(PoseidonHasher<E, S, R>);

impl<E: Engine, const S: usize, const R: usize> SpongePermutation<E> for PoseidonHasher<E, S, R> {
    fn permutation(&mut self) {
        debug_assert!(self.params.full_rounds & 1 == 0);
        let half_of_full_rounds = self.params.full_rounds / 2;

        let mut mds_result = [E::Fr::zero(); S];
        // full rounds
        for round in 0..half_of_full_rounds {
            // add round constatnts
            for (s, c) in self
                .state
                .iter_mut()
                .zip(&self.optimized_round_constants[round])
            {
                s.add_assign(c);
            }
            // apply sbox
            sbox::<E>(self.alpha, &mut self.state);
            // mul state by mds
            mmul_assign::<E, S>(&self.params.mds_matrix, &mut self.state);
        }

        // partial rounds
        // in this optimized version;
        // - first, use M' instead of sbox and matrix multiplication for other elements of state(not first element)
        // - second, instead of multiplication by original MDS matrix, multiply by M"(M" is a sparse matrix form)
        self.state
            .iter_mut()
            .zip(self.optimized_round_constants[half_of_full_rounds].iter())
            .for_each(|(s, c)| s.add_assign(c));
        mmul_assign::<E, S>(&self.optimized_mds_matrixes.0, &mut self.state);

        // this is an unrolled version of partial rounds
        for (round_constants, sparse_matrix) in self.optimized_round_constants
            [half_of_full_rounds + 1..half_of_full_rounds + self.params.partial_rounds]
            .iter()
            .chain(&[[E::Fr::zero(); S]])
            .zip(self.optimized_mds_matrixes.1.iter())
        {
            let mut quad = self.state[0];
            quad.square();
            quad.square();
            self.state[0].mul_assign(&quad);

            self.state[0].add_assign(&round_constants[0]);

            mds_result[0] = E::Fr::zero();
            for (a, b) in self.state.iter().zip(sparse_matrix[0].iter()) {
                let mut tmp = a.clone();
                tmp.mul_assign(&b);
                mds_result[0].add_assign(&tmp);
            }

            let mut tmp = sparse_matrix[1][0];
            tmp.mul_assign(&self.state[0]);
            tmp.add_assign(&self.state[1]);
            mds_result[1] = tmp;

            let mut tmp = sparse_matrix[2][0];
            tmp.mul_assign(&self.state[0]);
            tmp.add_assign(&self.state[2]);
            mds_result[2] = tmp;

            self.state.copy_from_slice(&mds_result[..]);
        }

        // full rounds
        for round in (self.params.partial_rounds + half_of_full_rounds)
            ..(self.params.partial_rounds + self.params.full_rounds)
        {
            // add round constants
            for (s, c) in self
                .state
                .iter_mut()
                .zip(&self.optimized_round_constants[round])
            {
                s.add_assign(c);
            }
            // apply sbox
            sbox::<E>(self.alpha, &mut self.state);

            // mul state by mds
            mmul_assign::<E, S>(&self.params.mds_matrix, &mut self.state);
        }
    }
}
