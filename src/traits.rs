
use franklin_crypto::{bellman::Engine};
#[derive(Debug, PartialEq, Eq)]
pub enum HashFamily {
    Rescue,
    Poseidon,
    RescuePrime,
}

pub trait HashParams<E: Engine, const RATE: usize, const WIDTH: usize>: Clone + Send + Sync {
    fn hash_family(&self) -> HashFamily;
    fn constants_of_round(&self, round: usize) -> [E::Fr; WIDTH];
    fn mds_matrix(&self) -> [[E::Fr; WIDTH]; WIDTH];
    fn number_of_full_rounds(&self) -> usize;
    fn number_of_partial_rounds(&self) -> usize;
    fn alpha(&self) -> E::Fr;
    fn alpha_inv(&self) -> E::Fr;
    fn optimized_round_constants(&self) -> &[[E::Fr; WIDTH]];
    fn optimized_mds_matrixes(&self) -> (&[[E::Fr; WIDTH]; WIDTH], &[[[E::Fr; WIDTH];WIDTH]]);
    fn can_use_custom_gates(&self) -> bool;
    fn set_allow_custom_gate(&mut self, allow: bool);
}