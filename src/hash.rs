use crate::traits::Sponge;
use crate::{common::domain_strategy::DomainStrategy, sponge::GenericSponge, traits::HashParams};
use franklin_crypto::bellman::Engine;
use std::convert::TryInto;

pub fn generic_hash<
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
    const LENGTH: usize,
>(
    params: &P,
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    inner_generic_hash(params, input, DomainStrategy::<RATE>::CustomFixedLength)
}

pub fn generic_hash_var_length<
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    params: &P,
    input: &[E::Fr],
) -> [E::Fr; RATE] {
    inner_generic_hash(params, input, DomainStrategy::<RATE>::CustomVariableLength)
}

fn inner_generic_hash<
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    params: &P,
    input: &[E::Fr],
    domain_strategy: DomainStrategy<RATE>,
) -> [E::Fr; RATE] {
    assert!(input.is_empty() == false, "empty input");
    let capacity_value = domain_strategy.compute_capacity::<E>(input.len());

    let padding_values = domain_strategy.generate_padding_values::<E>(input.len());

    let input_with_padding = input
        .iter()
        .chain(&padding_values[..])
        .cloned()
        .collect::<Vec<E::Fr>>();

    let mut sponge = GenericSponge::from_params(params);

    sponge.specialize(capacity_value);
    sponge.absorb(&input_with_padding);
    sponge.squeeze(None).try_into().expect("constant array")
}
