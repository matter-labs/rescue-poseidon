use crate::common::matrix::mmul_assign;
use crate::common::padding::PaddingStrategy;
use crate::sponge_impl;
use crate::{common::hash::generic_hash_with_padding, HasherParams};
use crate::{
    common::sbox::sbox,
    sponge::{SpongeParams, SpongePermutation, SpongeState, SpongeMode, SpongeModes, StatefulSponge},
};
use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;
/// This hash function received fixed length inputs and outputs
/// number of elements as equal to rate parameter.
/// This function apply fixed length padding strategy.
pub fn poseidon_hash_fixed_length<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, PoseidonHasher<E>>(
        input,
        PaddingStrategy::FixedLength(input.len()),
    )
}

/// A poseidon hash funciton outputs number of elements as equal to rate parameter
/// This function apply fixed length padding strategy.
pub fn poseidon_var_length<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, PoseidonHasher<E>>(
        input,
        PaddingStrategy::VariableLength(input.len()),
    )
}

/// This is hasher with a custom strategy which basically sets value of capacity
pub fn poseidon_hash<E: Engine>(input: &[E::Fr]) -> Vec<E::Fr> {
    generic_hash_with_padding::<E, PoseidonHasher<E>>(input, PaddingStrategy::Custom(input.len()))
}

#[derive(Clone, Debug)]
pub struct PoseidonHasher<E: Engine> {
    params: HasherParams<E>,
    state: Vec<E::Fr>,
    tmp_storage: Vec<E::Fr>,
    alpha: E::Fr,
    optimized_round_constants: Vec<Vec<E::Fr>>,
    optimized_mds_matrixes: (Vec<Vec<E::Fr>>, Vec<Vec<Vec<E::Fr>>>),
    sponge_mode: SpongeModes<E>,
}

impl<E: Engine> Default for PoseidonHasher<E> {
    fn default() -> Self {
        let (params, alpha, optimized_round_constants, optimized_mds_matrixes) =
            super::params::poseidon_light_params();
        Self {
            state: vec![E::Fr::zero(); params.state_width],
            tmp_storage: Vec::with_capacity(params.rate),
            alpha,
            params,
            optimized_round_constants,
            optimized_mds_matrixes,
            sponge_mode: SpongeModes::Standard,
        }
    }
}
impl<E: Engine> PoseidonHasher<E> {
    pub fn new(sponge_mode: SpongeModes<E>) -> Self {
        let (params, alpha, optimized_round_constants, optimized_mds_matrixes) =
            super::params::poseidon_light_params();
        Self {
            state: vec![E::Fr::zero(); params.state_width],
            tmp_storage: Vec::with_capacity(params.rate),
            alpha,
            params,
            optimized_round_constants,
            optimized_mds_matrixes,
            sponge_mode,
        }
    }
}

// Implementation of common parts
sponge_impl!(PoseidonHasher<E>);

impl<E: Engine> SpongePermutation<E> for PoseidonHasher<E> {
    fn permutation(&mut self) {
        debug_assert!(self.params.full_rounds & 1 == 0);
        let half_of_full_rounds = self.params.full_rounds / 2;

        let mut mds_result = vec![E::Fr::zero(); self.params.state_width];
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
            mmul_assign::<E>(&self.params.mds_matrix, &mut self.state);
        }

        // partial rounds
        // in this optimized version;
        // - first, use M' instead of sbox and matrix multiplication for other elements of state(not first element)
        // - second, instead of multiplication by original MDS matrix, multiply by M"(M" is a sparse matrix form)
        self.state
            .iter_mut()
            .zip(self.optimized_round_constants[half_of_full_rounds].iter())
            .for_each(|(s, c)| s.add_assign(c));
        mmul_assign::<E>(&self.optimized_mds_matrixes.0, &mut self.state);

        // this is an unrolled version of partial rounds
        for (round_constants, sparse_matrix) in self.optimized_round_constants
            [half_of_full_rounds + 1..half_of_full_rounds + self.params.partial_rounds]
            .iter()
            .chain(&[vec![E::Fr::zero(); self.params.state_width]])
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
            mmul_assign::<E>(&self.params.mds_matrix, &mut self.state);
        }
    }
}
