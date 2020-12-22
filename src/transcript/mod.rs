pub(crate) mod derived;

use crate::sponge::StatefulSponge;
use franklin_crypto::bellman::plonk::commitments::transcript::{Prng, Transcript};
use franklin_crypto::bellman::{Engine, PrimeField, PrimeFieldRepr};
use crate::stateful_transcript;

#[derive(Clone)]
pub struct RescueTranscript<E: Engine> {
    sponge: crate::rescue::RescueHasher<E>,
}
stateful_transcript!(RescueTranscript<E>, crate::rescue::RescueHasher::default);

#[derive(Clone)]
pub struct RescuePrimeTranscript<E: Engine> {
    sponge: crate::rescue_prime::RescuePrimeHasher<E>,
}
stateful_transcript!(RescuePrimeTranscript<E>, crate::rescue_prime::RescuePrimeHasher::default);

#[derive(Clone)]
pub struct PoseidonTranscript<E: Engine> {
    sponge: crate::poseidon::PoseidonHasher<E>,
}
stateful_transcript!(PoseidonTranscript<E>, crate::poseidon::PoseidonHasher::default);