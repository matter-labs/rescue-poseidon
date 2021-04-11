use super::sbox::*;
use super::utils::matrix_vector_product;
use crate::traits::{HashFamily, HashParams};
use franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;

use super::hash::{circuit_generic_hash, circuit_generic_hash_var_length};
use crate::rescue::RescueParams;
use franklin_crypto::bellman::SynthesisError;
use franklin_crypto::{
    bellman::Engine, plonk::circuit::allocated_num::Num,
    plonk::circuit::linear_combination::LinearCombination,
};

use std::convert::TryInto;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter.
/// Uses pre-defined state-width=3 and rate=2.
pub fn circuit_rescue_hash<E: Engine, CS: ConstraintSystem<E>, const L: usize>(
    cs: &mut CS,
    input: &[Num<E>; L],
) -> Result<[Num<E>; 2], SynthesisError> {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let params = RescueParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash(cs, &params, input).map(|res| res.try_into().expect(""))
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule.
/// Uses pre-defined state-width=3 and rate=2.
pub fn gadget_rescue_hash_var_length<E: Engine, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<[Num<E>; 2], SynthesisError> {
    // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let params = RescueParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash_var_length(cs, &params, input).map(|res| res.try_into().expect(""))
}

pub fn gadget_generic_rescue_hash<
    E: Engine,
    CS: ConstraintSystem<E>,
    const RATE: usize,
    const WIDTH: usize,
    const LENGTH: usize,
>(
    cs: &mut CS,
    input: &[Num<E>; LENGTH],
) -> Result<[Num<E>; RATE], SynthesisError> {
    let params = RescueParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash(cs, &params, input).map(|res| res.try_into().expect(""))
}

pub fn gadget_generic_rescue_hash_var_length<
    E: Engine,
    CS: ConstraintSystem<E>,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<[Num<E>; RATE], SynthesisError> {
    let params = RescueParams::<E, RATE, WIDTH>::default();
    circuit_generic_hash_var_length(cs, &params, input).map(|res| res.try_into().expect(""))
}

pub(crate) fn gadget_rescue_round_function<
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
        HashFamily::Rescue,
        "Incorrect hash family!"
    );
    state
        .iter_mut()
        .zip(params.constants_of_round(0).iter())
        .for_each(|(s, c)| s.add_assign_constant(*c));

    for round in 0..2 * params.number_of_full_rounds() {
        // apply sbox
        if round & 1 == 0 {
            sbox_quintic_inv::<E, _>(cs, params.alpha_inv(), state)?;
        } else {
            sbox_quintic(cs, state)?;
        }
        // mds row
        // TODO remove mut from mds
        *state = matrix_vector_product(cs, &params.mds_matrix(), state)?;

        // round constants
        for (s, c) in state
            .iter_mut()
            .zip(params.constants_of_round(round + 1).iter().cloned())
        {
            s.add_assign_constant(c);
        }
    }
    Ok(())
}
