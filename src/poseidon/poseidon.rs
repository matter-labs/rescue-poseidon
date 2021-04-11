use crate::common::{matrix::mmul_assign, sbox::sbox};
use crate::hash::{generic_hash, generic_hash_var_length};
use crate::traits::{HashFamily, HashParams};
use franklin_crypto::bellman::{Engine, Field};

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter. Uses state-width=3 and rate=2.
pub fn poseidon_hash<E: Engine, const L: usize>(input: &[E::Fr; L]) -> [E::Fr; 2] {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let params = PoseidonParams::<E, RATE, WIDTH>::default();
    generic_hash(&params, input)
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule. Uses state-width=3 and rate=2.
pub fn poseidon_hash_var_length<E: Engine>(input: &[E::Fr]) -> [E::Fr; 2] {
    // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let params = PoseidonParams::<E, RATE, WIDTH>::default();
    generic_hash_var_length(&params, input)
}

pub fn generic_poseidon_hash<
    E: Engine,
    const RATE: usize,
    const WIDTH: usize,
    const LENGTH: usize,
>(
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    let params = PoseidonParams::<E, RATE, WIDTH>::default();
    generic_hash(&params, input)
}

pub fn generic_poseidon_hash_var_length<
    E: Engine,
    const RATE: usize,
    const WIDTH: usize,
>(
    input: &[E::Fr],
) -> [E::Fr; RATE] {
    let params = PoseidonParams::<E, RATE, WIDTH>::default();
    generic_hash_var_length(&params, input)
}

#[derive(Clone, Debug)]
pub struct PoseidonParams<E: Engine, const RATE: usize, const WIDTH: usize> {
    state: [E::Fr; WIDTH],
    mds_matrix: [[E::Fr; WIDTH]; WIDTH],
    optimized_round_constants: Vec<[E::Fr; WIDTH]>,
    optimized_mds_matrixes: (
        [[E::Fr; WIDTH]; WIDTH],
        Vec<[[E::Fr; WIDTH]; WIDTH]>,
    ),
    alpha: E::Fr,
    full_rounds: usize,
    partial_rounds: usize,
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> Default
    for PoseidonParams<E, RATE, WIDTH>
{
    fn default() -> Self {
        let (params, alpha, optimized_round_constants, optimized_mds_matrixes) =
            super::params::poseidon_light_params::<E, RATE, WIDTH>();
        Self {
            state: [E::Fr::zero(); WIDTH],
            mds_matrix: params.mds_matrix,
            alpha,
            optimized_round_constants,
            optimized_mds_matrixes,
            full_rounds: params.full_rounds,
            partial_rounds: params.partial_rounds,
        }
    }
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> HashParams<E, RATE, WIDTH>
    for PoseidonParams<E, RATE, WIDTH>
{
    fn hash_family(&self) -> HashFamily {
        HashFamily::Poseidon
    }

    fn constants_of_round(&self, _round: usize) -> [E::Fr; WIDTH] {
        unimplemented!("Poseidon uses optimized constants")
    }

    fn mds_matrix(&self) -> [[E::Fr; WIDTH]; WIDTH] {
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

    fn optimized_round_constants(&self) -> &[[E::Fr; WIDTH]] {
        &self.optimized_round_constants
    }

    fn optimized_mds_matrixes(
        &self,
    ) -> (
        &[[E::Fr; WIDTH]; WIDTH],
        &[[[E::Fr; WIDTH]; WIDTH]],
    ) {
        (
            &self.optimized_mds_matrixes.0,
            &self.optimized_mds_matrixes.1,
        )
    }
}

pub(crate) fn poseidon_round_function<
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    params: &P,
    state: &mut [E::Fr; WIDTH],
) {
    assert_eq!(params.hash_family(), HashFamily::Poseidon, "Incorrect hash family!");
    debug_assert!(params.number_of_full_rounds() & 1 == 0);
    let half_of_full_rounds = params.number_of_full_rounds() / 2;

    let mut mds_result = [E::Fr::zero(); WIDTH];

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
        mmul_assign::<E, WIDTH>(&params.mds_matrix(), state);
    }

    // partial rounds
    // in this optimized version;
    // - first, use M' instead of sbox and matrix multiplication for other elements of state(not first element)
    // - second, instead of multiplication by original MDS matrix, multiply by M"(M" is a sparse matrix form)

    state
        .iter_mut()
        .zip(optimized_round_constants[half_of_full_rounds].iter())
        .for_each(|(s, c)| s.add_assign(c));
    mmul_assign::<E, WIDTH>(&sparse_matrixes.0, state);

    // this is an unrolled version of partial rounds
    for (round_constants, sparse_matrix) in optimized_round_constants
        [half_of_full_rounds + 1..half_of_full_rounds + params.number_of_partial_rounds()]
        .iter()
        .chain(&[[E::Fr::zero(); WIDTH]])
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
        mmul_assign::<E, WIDTH>(&params.mds_matrix(), state);
    }
}
