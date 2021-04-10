use crate::common::{matrix::mmul_assign, sbox::sbox};
use crate::hash::{generic_hash, generic_hash_var_length};
use crate::traits::{HashFamily, HashParams};
use franklin_crypto::bellman::{Engine, Field};

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter. Uses state-width=3 and rate=2.
pub fn poseidon_hash<E: Engine, const L: usize>(input: &[E::Fr; L]) -> [E::Fr; 2] {
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;

    let params = PoseidonParams::<E, STATE_WIDTH, RATE>::default();
    generic_hash(&params, input)
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule. Uses state-width=3 and rate=2.
pub fn poseidon_hash_var_length<E: Engine>(input: &[E::Fr]) -> [E::Fr; 2] {
    // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;

    let params = PoseidonParams::<E, STATE_WIDTH, RATE>::default();
    generic_hash_var_length(&params, input)
}

pub fn generic_poseidon_hash<
    E: Engine,
    const STATE_WIDTH: usize,
    const RATE: usize,
    const LENGTH: usize,
>(
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    let params = PoseidonParams::<E, STATE_WIDTH, RATE>::default();
    generic_hash(&params, input)
}

pub fn generic_poseidon_hash_var_length<
    E: Engine,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    input: &[E::Fr],
) -> [E::Fr; RATE] {
    let params = PoseidonParams::<E, STATE_WIDTH, RATE>::default();
    generic_hash_var_length(&params, input)
}

#[derive(Clone, Debug)]
pub struct PoseidonParams<E: Engine, const STATE_WIDTH: usize, const RATE: usize> {
    state: [E::Fr; STATE_WIDTH],
    mds_matrix: [[E::Fr; STATE_WIDTH]; STATE_WIDTH],
    optimized_round_constants: Vec<[E::Fr; STATE_WIDTH]>,
    optimized_mds_matrixes: (
        [[E::Fr; STATE_WIDTH]; STATE_WIDTH],
        Vec<[[E::Fr; STATE_WIDTH]; STATE_WIDTH]>,
    ),
    alpha: E::Fr,
    full_rounds: usize,
    partial_rounds: usize,
}

impl<E: Engine, const STATE_WIDTH: usize, const RATE: usize> Default
    for PoseidonParams<E, STATE_WIDTH, RATE>
{
    fn default() -> Self {
        let (params, alpha, optimized_round_constants, optimized_mds_matrixes) =
            super::params::poseidon_light_params::<E, STATE_WIDTH, RATE>();
        Self {
            state: [E::Fr::zero(); STATE_WIDTH],
            mds_matrix: params.mds_matrix,
            alpha,
            optimized_round_constants,
            optimized_mds_matrixes,
            full_rounds: params.full_rounds,
            partial_rounds: params.partial_rounds,
        }
    }
}

impl<E: Engine, const STATE_WIDTH: usize, const RATE: usize> HashParams<E, STATE_WIDTH, RATE>
    for PoseidonParams<E, STATE_WIDTH, RATE>
{
    fn hash_family(&self) -> HashFamily {
        HashFamily::Poseidon
    }

    fn constants_of_round(&self, _round: usize) -> [E::Fr; STATE_WIDTH] {
        unimplemented!("Poseidon uses optimized constants")
    }

    fn mds_matrix(&self) -> [[E::Fr; STATE_WIDTH]; STATE_WIDTH] {
        self.mds_matrix
    }

    fn number_of_full_rounds(&self) -> usize {
        self.full_rounds
    }

    fn number_of_partial_rounds(&self) -> usize {
        self.partial_rounds
    }

    fn alpha(&self) -> E::Fr {
        self.alpha
    }

    fn alpha_inv(&self) -> E::Fr {
        unimplemented!("Poseidon doesn't have inverse direction")
    }

    fn optimized_round_constants(&self) -> &[[E::Fr; STATE_WIDTH]] {
        &self.optimized_round_constants
    }

    fn optimized_mds_matrixes(
        &self,
    ) -> (
        &[[E::Fr; STATE_WIDTH]; STATE_WIDTH],
        &[[[E::Fr; STATE_WIDTH]; STATE_WIDTH]],
    ) {
        (
            &self.optimized_mds_matrixes.0,
            &self.optimized_mds_matrixes.1,
        )
    }
}

pub(crate) fn poseidon_round_function<
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    params: &P,
    state: &mut [E::Fr; STATE_WIDTH],
) {
    assert_eq!(params.hash_family(), HashFamily::Poseidon, "Incorrect hash family!");
    debug_assert!(params.number_of_full_rounds() & 1 == 0);
    let half_of_full_rounds = params.number_of_full_rounds() / 2;

    let mut mds_result = [E::Fr::zero(); STATE_WIDTH];

    let optimized_round_constants = params.optimized_round_constants();
    let sparse_matrixes = params.optimized_mds_matrixes();
    // full rounds
    for round in 0..half_of_full_rounds {
        // add round constatnts
        for (s, c) in state.iter_mut().zip(&optimized_round_constants[round]) {
            s.add_assign(c);
        }
        // apply sbox
        sbox::<E>(params.alpha(), state);
        // mul state by mds
        mmul_assign::<E, STATE_WIDTH>(&params.mds_matrix(), state);
    }

    // partial rounds
    // in this optimized version;
    // - first, use M' instead of sbox and matrix multiplication for other elements of state(not first element)
    // - second, instead of multiplication by original MDS matrix, multiply by M"(M" is a sparse matrix form)

    state
        .iter_mut()
        .zip(optimized_round_constants[half_of_full_rounds].iter())
        .for_each(|(s, c)| s.add_assign(c));
    mmul_assign::<E, STATE_WIDTH>(&sparse_matrixes.0, state);

    // this is an unrolled version of partial rounds
    for (round_constants, sparse_matrix) in optimized_round_constants
        [half_of_full_rounds + 1..half_of_full_rounds + params.number_of_partial_rounds()]
        .iter()
        .chain(&[[E::Fr::zero(); STATE_WIDTH]])
        .zip(sparse_matrixes.1.iter())
    {
        let mut quad = state[0];
        quad.square();
        quad.square();
        state[0].mul_assign(&quad);

        state[0].add_assign(&round_constants[0]);

        mds_result[0] = E::Fr::zero();
        for (a, b) in state.iter().zip(sparse_matrix[0].iter()) {
            let mut tmp = a.clone();
            tmp.mul_assign(&b);
            mds_result[0].add_assign(&tmp);
        }

        let mut tmp = sparse_matrix[1][0];
        tmp.mul_assign(&state[0]);
        tmp.add_assign(&state[1]);
        mds_result[1] = tmp;

        let mut tmp = sparse_matrix[2][0];
        tmp.mul_assign(&state[0]);
        tmp.add_assign(&state[2]);
        mds_result[2] = tmp;

        state.copy_from_slice(&mds_result[..]);
    }

    // full rounds
    for round in (params.number_of_partial_rounds() + half_of_full_rounds)
        ..(params.number_of_partial_rounds() + params.number_of_full_rounds())
    {
        // add round constants
        for (s, c) in state.iter_mut().zip(&optimized_round_constants[round]) {
            s.add_assign(c);
        }
        // apply sbox
        sbox::<E>(params.alpha(), state);

        // mul state by mds
        mmul_assign::<E, STATE_WIDTH>(&params.mds_matrix(), state);
    }
}
