use super::domain_strategy::DomainStrategy;
use crate::sponge::StatefulSponge;
use franklin_crypto::bellman::Engine;

// This function specializes sponge construction with respect to domain strategy.
pub(crate) fn generic_hash_with_padding<E: Engine, Sponge: StatefulSponge<E, S, R>, const S: usize, const R: usize>(
    input: &[E::Fr],
    padding_strategy: DomainStrategy<R>,
) -> Vec<E::Fr> {
    let mut sponge = Sponge::default();
    
    let capacity_value = padding_strategy.compute_capacity::<E>(input.len());

    let padding_values = padding_strategy.generate_padding_values::<E>(input.len());

    let input_with_padding = input
        .iter()
        .chain(&padding_values[..])
        .cloned()
        .collect::<Vec<E::Fr>>();

    sponge.specialize(capacity_value);
    sponge.absorb(&input_with_padding);
    sponge.squeeze(None)
}
