use crate::common::matrix::mmul_assign;
use crate::common::sbox::sbox;
use crate::hash::{generic_hash, generic_hash_var_length};
use crate::traits::{HashFamily, HashParams};
use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;
use std::convert::TryInto;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter.
/// Uses pre-defined state-width=3 and rate=2.
pub fn rescue_prime_hash<E: Engine, const L: usize>(input: &[E::Fr; L]) -> [E::Fr; 2] {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    generic_hash(&params, input)
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule.
/// Uses pre-defined state-width=3 and rate=2.
pub fn rescue_prime_hash_var_length<E: Engine>(input: &[E::Fr]) -> [E::Fr; 2] {
    // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    generic_hash_var_length(&params, input)
}

pub fn generic_rescue_prime<
    E: Engine,
    const RATE: usize,
    const WIDTH: usize,
    const LENGTH: usize,
>(
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    generic_hash(&params, input)
}

pub fn generic_rescue_prime_var_length<E: Engine, const RATE: usize, const WIDTH: usize>(
    input: &[E::Fr],
) -> [E::Fr; RATE] {
    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    generic_hash_var_length(&params, input)
}

#[derive(Clone, Debug)]
pub struct RescuePrimeParams<E: Engine, const RATE: usize, const WIDTH: usize> {
    pub full_rounds: usize,
    pub round_constants: Vec<[E::Fr; WIDTH]>,
    pub mds_matrix: [[E::Fr; WIDTH]; WIDTH],
    pub alpha: E::Fr,
    pub alpha_inv: E::Fr,
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> Default
    for RescuePrimeParams<E, RATE, WIDTH>
{
    fn default() -> Self {
        let (params, alpha, alpha_inv) = super::params::rescue_prime_params::<E, RATE, WIDTH>();
        Self {
            full_rounds: params.full_rounds,
            round_constants: params.round_constants().try_into().expect("constant array"),
            mds_matrix: *params.mds_matrix(),
            alpha,
            alpha_inv,
        }
    }
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> HashParams<E, RATE, WIDTH>
    for RescuePrimeParams<E, RATE, WIDTH>
{
    fn hash_family(&self) -> HashFamily {
        HashFamily::RescuePrime
    }

    fn constants_of_round(&self, round: usize) -> [E::Fr; WIDTH] {
        self.round_constants[round]
    }

    fn mds_matrix(&self) -> [[E::Fr; WIDTH]; WIDTH] {
        self.mds_matrix
    }

    fn number_of_full_rounds(&self) -> usize {
        self.full_rounds
    }

    fn number_of_partial_rounds(&self) -> usize {
        unimplemented!("RescuePrime doesn't have partial rounds.")
    }

    fn alpha(&self) -> E::Fr {
        self.alpha
    }

    fn alpha_inv(&self) -> E::Fr {
        self.alpha_inv
    }

    fn optimized_mds_matrixes(&self) -> (&[[E::Fr; WIDTH]; WIDTH], &[[[E::Fr; WIDTH]; WIDTH]]) {
        unimplemented!("RescuePrime doesn't use optimized mds matrixes")
    }

    fn optimized_round_constants(&self) -> &[[E::Fr; WIDTH]] {
        unimplemented!("RescuePrime doesn't use optimized round constants")
    }
}

pub(crate) fn rescue_prime_round_function<
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    params: &P,
    state: &mut [E::Fr; WIDTH],
) {
    assert_eq!(
        params.hash_family(),
        HashFamily::RescuePrime,
        "Incorrect hash family!"
    );
    for round in 0..params.number_of_full_rounds() - 1 {
        // sbox alpha
        sbox::<E>(params.alpha(), state);
        // mds
        mmul_assign::<E, WIDTH>(&params.mds_matrix(), state);

        // round constants
        state
            .iter_mut()
            .zip(params.constants_of_round(round).iter())
            .for_each(|(s, c)| s.add_assign(c));
        // sbox alpha inv
        sbox::<E>(params.alpha_inv(), state);

        // mds
        mmul_assign::<E, WIDTH>(&params.mds_matrix(), state);

        // round constants
        state
            .iter_mut()
            .zip(params.constants_of_round(round + 1).iter())
            .for_each(|(s, c)| s.add_assign(c));
    }
}
