use crate::common::params::HasherParams;
use crate::common::{
    domain_strategy::DomainStrategy, hash::generic_hash_with_padding, matrix::mmul_assign,
    sbox::sbox,
};
use crate::sponge::{SpongeMode, SpongeModes, SpongePermutation, SpongeState, StatefulSponge};
use crate::sponge_impl;
use franklin_crypto::bellman::{Engine, Field};
use std::convert::TryInto;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter. Uses state-width=3 and rate=2.
pub fn rescue_hash<E: Engine, const L: usize>(input: &[E::Fr; L]) -> [E::Fr; 2] {
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;

    rescue_generic_fixed_length::<E, STATE_WIDTH, RATE, L>(input)
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule. Uses state-width=3 and rate=2.
pub fn rescue_hash_var_length<E: Engine>(input: &[E::Fr]) -> [E::Fr; 2] {
    // TODO: try to implement const_generics_defaults: https://github.com/rust-lang/rust/issues/44580
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;

    rescue_generic_var_length::<E, STATE_WIDTH, RATE>(input)
}

pub(crate) fn rescue_generic_fixed_length<
    E: Engine,
    const STATE_WIDTH: usize,
    const RATE: usize,
    const LENGTH: usize,
>(
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    let result =
        generic_hash_with_padding::<E, RescueHasher<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            input,
            DomainStrategy::CustomFixedLength,
        );

    result.try_into().expect("fixed length array")
}

pub(crate) fn rescue_generic_var_length<E: Engine, const STATE_WIDTH: usize, const RATE: usize>(
    input: &[E::Fr],
) -> [E::Fr; RATE] {
    let result =
        generic_hash_with_padding::<E, RescueHasher<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            input,
            DomainStrategy::CustomVariableLength,
        );

    result.try_into().expect("fixed length array")
}

pub trait HashParams<E: Engine, const STATE_WIDTH: usize, const RATE: usize> {
    fn constants_of_round(&self, round: usize) -> [E::Fr; STATE_WIDTH];
    fn mds_matrix(&self) -> [[E::Fr; STATE_WIDTH]; STATE_WIDTH];
    fn number_of_full_rounds(&self) -> usize;
    fn number_of_partial_rounds(&self) -> usize;
    fn alpha(&self) -> E::Fr;
    fn alpha_inv(&self) -> E::Fr;
}
#[derive(Clone, Debug)]
pub struct RescueParams<E: Engine, const STATE_WIDTH: usize, const RATE: usize> {
    pub full_rounds: usize,
    pub round_constants: Vec<[E::Fr; STATE_WIDTH]>,
    pub mds_matrix: [[E::Fr; STATE_WIDTH]; STATE_WIDTH],
    pub alpha: E::Fr,
    pub alpha_inv: E::Fr,
}

impl<E: Engine, const STATE_WIDTH: usize, const RATE: usize> Default
    for RescueParams<E, STATE_WIDTH, RATE>
{
    fn default() -> Self {
        let (params, alpha, alpha_inv) = super::params::rescue_params::<E, STATE_WIDTH, RATE>();
        Self {
            full_rounds: params.full_rounds,
            round_constants: params
                .round_constants()
                .try_into()
                .expect("round constants"),
            mds_matrix: *params.mds_matrix(),
            alpha,
            alpha_inv,
        }
    }
}

impl<E: Engine, const STATE_WIDTH: usize, const RATE: usize> HashParams<E, STATE_WIDTH, RATE>
    for RescueParams<E, STATE_WIDTH, RATE>
{
    fn constants_of_round(&self, round: usize) -> [E::Fr; STATE_WIDTH] {
        self.round_constants[round]
    }

    fn mds_matrix(&self) -> [[E::Fr; STATE_WIDTH]; STATE_WIDTH] {
        self.mds_matrix
    }

    fn number_of_full_rounds(&self) -> usize {
        self.full_rounds
    }

    fn number_of_partial_rounds(&self) -> usize {
        unimplemented!("Rescue doesn't have partial rounds.")
    }

    fn alpha(&self) -> E::Fr {
        self.alpha
    }

    fn alpha_inv(&self) -> E::Fr {
        self.alpha_inv
    }
}

impl<E: Engine, const STATE_WIDTH: usize, const RATE: usize> RescueParams<E, STATE_WIDTH, RATE> {
    pub fn get_hasher() -> RescueHasher<E, STATE_WIDTH, RATE> {
        RescueHasher::new_from_params(RescueParams::default())
    }
}

#[derive(Debug, Clone)]
pub struct RescueHasher<E: Engine, const STATE_WIDTH: usize, const RATE: usize> {
    params: RescueParams<E, STATE_WIDTH, RATE>,
    state: [E::Fr; STATE_WIDTH],
    sponge_mode: SpongeModes,
}

impl<E: Engine, const STATE_WIDTH: usize, const RATE: usize> Default
    for RescueHasher<E, STATE_WIDTH, RATE>
{
    fn default() -> Self {
        Self {
            state: [E::Fr::zero(); STATE_WIDTH],
            sponge_mode: SpongeModes::Standard(false),
            params: RescueParams::default(),
        }
    }
}

impl<E: Engine, const STATE_WIDTH: usize, const RATE: usize> RescueHasher<E, STATE_WIDTH, RATE> {
    pub fn new_from_params(params: RescueParams<E, STATE_WIDTH, RATE>) -> Self {
        Self {
            state: [E::Fr::zero(); STATE_WIDTH],
            sponge_mode: SpongeModes::Duplex(false),
            params,
        }
    }
    pub fn new_duplex() -> Self {
        Self {
            state: [E::Fr::zero(); STATE_WIDTH],
            sponge_mode: SpongeModes::Duplex(false),
            params: RescueParams::default(),
        }
    }
}

// common parts of sponge
sponge_impl!(RescueHasher<E, S, R>);

impl<E: Engine, const S: usize, const R: usize> SpongePermutation<E> for RescueHasher<E, S, R> {
    fn permutation(&mut self) {
        rescue_round_function(&self.params, &mut self.state)
    }
}

pub fn rescue_round_function<
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    params: &P,
    state: &mut [E::Fr; STATE_WIDTH],
) {
    // round constants for first step
    state
        .iter_mut()
        .zip(params.constants_of_round(0).iter())
        .for_each(|(s, c)| s.add_assign(c));

    for round in 0..2 * params.number_of_full_rounds() {
        // sbox
        if round & 1 == 0 {
            sbox::<E>(params.alpha_inv(), state);
        } else {
            sbox::<E>(params.alpha(), state);
        }

        // mds
        mmul_assign::<E, STATE_WIDTH>(&params.mds_matrix(), state);

        // round constants
        state
            .iter_mut()
            .zip(params.constants_of_round(round + 1).iter())
            .for_each(|(s, c)| s.add_assign(c));
    }
}
