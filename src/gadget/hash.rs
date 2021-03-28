use super::sponge::StatefulSpongeGadget;
use crate::common::padding::PaddingStrategy;
use franklin_crypto::{bellman::{Engine, SynthesisError}, plonk::circuit::linear_combination::LinearCombination};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::allocated_num::Num,
};

// This function specializes sponge construction, computes required
// values for padding and outputs hash of pre-image. 
pub(crate) fn generic_hash<E: Engine, CS: ConstraintSystem<E>, SPONGE: StatefulSpongeGadget<E, R, S>, const R: usize, const S: usize>(
    cs: &mut CS,
    input: &[Num<E>],
    padding_strategy: PaddingStrategy,
) -> Result<Vec<Num<E>>, SynthesisError> {
    let mut sponge = SPONGE::default();
    let rate = R;

    let capacity_value = padding_strategy.compute_capacity::<E>(input.len(), rate).map(|el| {
        let mut lc = LinearCombination::zero();
        lc.add_assign_constant(el);
        lc
    });
    let padding_values = padding_strategy
        .generate_padding_values::<E>(input.len(), rate)
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

    // TODO
    let output = sponge.squeeze(cs, None)?;
    Ok(output)
}