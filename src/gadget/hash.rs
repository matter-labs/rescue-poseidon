use super::sponge::StatefulSpongeGadget;
use crate::common::domain_strategy::DomainStrategy;
use franklin_crypto::{bellman::{Engine, SynthesisError}, plonk::circuit::linear_combination::LinearCombination};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::allocated_num::Num,
};

// This function specializes sponge construction with respect to domain strategy.
pub(crate) fn generic_hash<E: Engine, CS: ConstraintSystem<E>, SPONGE: StatefulSpongeGadget<E, S, R>, const S: usize, const R: usize>(
    cs: &mut CS,
    input: &[Num<E>],
    domain_strategy: DomainStrategy<R>,
) -> Result<Vec<Num<E>>, SynthesisError> {
    let mut sponge = SPONGE::default();

    let capacity_value = domain_strategy.compute_capacity::<E>(input.len()).map(|el| {
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


    sponge.specialize(capacity_value);
    sponge.absorb(cs, &input_with_padding)?;

    let output = sponge.squeeze(cs, None)?;
    Ok(output)
}