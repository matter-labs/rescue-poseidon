use super::sponge::StatefulSpongeGadget;
use crate::common::padding::PaddingStrategy;
use franklin_crypto::{bellman::{Engine, SynthesisError}, plonk::circuit::linear_combination::LinearCombination};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::allocated_num::Num,
};

// This function specializes sponge construction, computes required
// values for padding and outputs hash of pre-image. 
pub(crate) fn generic_hash<E: Engine, CS: ConstraintSystem<E>, S: StatefulSpongeGadget<E>>(
    cs: &mut CS,
    input: &[Num<E>],
    padding_strategy: PaddingStrategy,
) -> Result<Vec<Num<E>>, SynthesisError> {
    let mut sponge = S::default();
    let rate = sponge.rate();

    let capacity_value = padding_strategy.compute_capacity::<E>(rate).map(|el| {
        let mut lc = LinearCombination::zero();
        lc.add_assign_constant(el);
        lc
    });
    let padding_values = padding_strategy
        .generate_padding_values::<E>(rate)
        .iter()
        .map(|el| Num::Constant(*el))
        .collect::<Vec<Num<E>>>();

    let input_with_padding = input
        .iter()
        .chain(&padding_values[..])
        .cloned()
        .collect::<Vec<Num<E>>>();

    sponge.absorb_multi(cs, &input)?;

    sponge.specialize(capacity_value);
    sponge.absorb_multi(cs, &input_with_padding)?;

    let output = sponge.squeeze(cs)?;
    Ok(output)
}