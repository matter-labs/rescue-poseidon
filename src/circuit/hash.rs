use super::sponge::GenericSpongeGadget;
use super::traits::SpongeGadget;
use crate::{common::domain_strategy::DomainStrategy, traits::HashParams};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::allocated_num::Num,
};
use franklin_crypto::{
    bellman::{Engine, SynthesisError},
    plonk::circuit::linear_combination::LinearCombination,
};
use std::convert::TryInto;

pub fn circuit_generic_hash<
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    CS: ConstraintSystem<E>,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    cs: &mut CS,
    params: &P,
    input: &[Num<E>],
) -> Result<Vec<Num<E>>, SynthesisError> {
    inner_generic_hash(cs, params, input, DomainStrategy::<RATE>::CustomFixedLength)
}

pub fn circuit_generic_hash_var_length<
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    CS: ConstraintSystem<E>,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    cs: &mut CS,
    params: &P,
    input: &[Num<E>],
) -> Result<Vec<Num<E>>, SynthesisError> {
    inner_generic_hash(cs, params, input, DomainStrategy::<RATE>::CustomVariableLength)
}


fn inner_generic_hash<
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    CS: ConstraintSystem<E>,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    cs: &mut CS,
    params: &P,
    input: &[Num<E>],
    domain_strategy: DomainStrategy<RATE>,
) -> Result<Vec<Num<E>>, SynthesisError> {
    let capacity_value = domain_strategy
        .compute_capacity::<E>(input.len())
        .map(|el| {
            let mut lc = LinearCombination::zero();
            lc.add_assign_constant(el);
            lc
        });

    let padding_values = domain_strategy
        .generate_padding_values::<E>(input.len())
        .iter()
        .map(|el| Num::Constant(*el))
        .collect::<Vec<Num<E>>>();

    let input_with_padding = input
        .iter()
        .chain(&padding_values[..])
        .cloned()
        .collect::<Vec<Num<E>>>();

    let mut sponge = GenericSpongeGadget::from(params);

    sponge.specialize(capacity_value)?;
    sponge.absorb(cs, &input_with_padding)?;
    let output = sponge
        .squeeze(cs, None)?
        .try_into()
        .expect("constant array");

    Ok(output)
}
