use super::padding::PaddingStrategy;
use crate::sponge::StatefulSponge;
use franklin_crypto::bellman::Engine;

// This function specializes sponge construction, computes required
// values for padding and outputs hash of pre-image. 
pub(crate) fn generic_hash_with_padding<E: Engine, Sponge: StatefulSponge<E, S, R>, const S: usize, const R: usize>(
    input: &[E::Fr],
    padding_strategy: PaddingStrategy,
) -> Vec<E::Fr> {
    let mut sponge = Sponge::default();
    let rate = R;
    
    let capacity_value = padding_strategy.compute_capacity::<E>(input.len(), rate);
    let padding_values = padding_strategy.generate_padding_values::<E>(input.len(), rate);

    let input_with_padding = input
        .iter()
        .chain(&padding_values[..])
        .cloned()
        .collect::<Vec<E::Fr>>();

    sponge.specialize(capacity_value);
    // sponge.absorb_multi(&input_with_padding);
    // TODO
    sponge.absorb(&input_with_padding);
    sponge.squeeze(None) // TODO: what is the length of output?
}
