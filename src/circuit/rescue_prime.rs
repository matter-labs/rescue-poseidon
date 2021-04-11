use super::hash::{circuit_generic_hash, circuit_generic_hash_var_length};
use super::sbox::*;
use super::utils::matrix_vector_product;
use crate::rescue_prime::RescuePrimeParams;
use crate::traits::{HashFamily, HashParams};
use franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use franklin_crypto::bellman::SynthesisError;
use franklin_crypto::{
    bellman::Engine,
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};
use std::convert::TryInto;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter.
/// Uses pre-defined state-width=3 and rate=2.
pub fn gadget_rescue_prime_hash<E: Engine, CS: ConstraintSystem<E>, const L: usize>(
    cs: &mut CS,
    input: &[Num<E>; L],
) -> Result<[Num<E>; 2], SynthesisError> {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash(cs, &params, input).map(|res| res.try_into().expect(""))
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule.
/// Uses pre-defined state-width=3 and rate=2.
pub fn gadget_rescue_prime_hash_var_length<E: Engine, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<[Num<E>; 2], SynthesisError> {
    // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash_var_length(cs, &params, input).map(|res| res.try_into().expect(""))
}

pub fn gadget_generic_rescue_prime_hash<
    E: Engine,
    CS: ConstraintSystem<E>,
    const RATE: usize,
    const WIDTH: usize,
    const LENGTH: usize,
>(
    cs: &mut CS,
    input: &[Num<E>; LENGTH],
) -> Result<[Num<E>; RATE], SynthesisError> {
    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash(cs, &params, input).map(|res| res.try_into().expect(""))
}

pub fn gadget_generic_rescue_prime_hash_var_length<
    E: Engine,
    CS: ConstraintSystem<E>,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<[Num<E>; RATE], SynthesisError> {
    let params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash_var_length(cs, &params, input).map(|res| res.try_into().expect(""))
}

pub(crate) fn gadget_rescue_prime_round_function<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    params: &P,
    state: &mut [LinearCombination<E>; WIDTH],
) -> Result<(), SynthesisError> {
    assert_eq!(
        params.hash_family(),
        HashFamily::RescuePrime,
        "Incorrect hash family!"
    );

    for round in 0..params.number_of_full_rounds() {
        // apply sbox
        // each lc will have 3 terms but there will be 1 in first iteration
        // total cost 2 gate per state vars = 6
        sbox_quintic(cs, state)?;

        // mul by mds
        *state = matrix_vector_product(cs, &params.mds_matrix(), state)?;

        // round constants
        let constants = params.constants_of_round(round);
        for (s, c) in state.iter_mut().zip(constants.iter().cloned()) {
            s.add_assign_constant(c);
        }
        // apply inverse sbox
        sbox_quintic_inv::<E, _>(cs, params.alpha_inv(), state)?;

        // mul by mds
        *state = matrix_vector_product(cs, &params.mds_matrix(), state)?;

        // round constants
        let constants = params.constants_of_round(round + 1);
        for (s, c) in state.iter_mut().zip(constants.iter().cloned()) {
            s.add_assign_constant(c);
        }
    }
    Ok(())
}
