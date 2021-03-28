pub(crate) mod derived;

use crate::sponge::StatefulSponge;
use crate::stateful_transcript;
use franklin_crypto::bellman::plonk::commitments::transcript::{Prng, Transcript};
use franklin_crypto::bellman::{Engine, PrimeField, PrimeFieldRepr};

#[derive(Clone)]
pub struct RescueTranscript<E: Engine, const S: usize, const R: usize> {
    sponge: crate::rescue::RescueHasher<E, S, R>,
}
stateful_transcript!(RescueTranscript<E,S,R>, crate::rescue::RescueHasher::default);

#[derive(Clone)]
pub struct RescuePrimeTranscript<E: Engine, const S: usize, const R: usize> {
    sponge: crate::rescue_prime::RescuePrimeHasher<E, S, R>,
}
stateful_transcript!(RescuePrimeTranscript<E,S,R>, crate::rescue_prime::RescuePrimeHasher::default);

#[derive(Clone)]
pub struct PoseidonTranscript<E: Engine, const S: usize, const R: usize> {
    sponge: crate::poseidon::PoseidonHasher<E, S, R>,
}
stateful_transcript!(PoseidonTranscript<E, S, R>, crate::poseidon::PoseidonHasher::default);
