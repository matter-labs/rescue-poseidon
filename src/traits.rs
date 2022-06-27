use franklin_crypto::bellman::Engine;
#[derive(Debug, PartialEq, Eq)]
pub enum HashFamily {
    Rescue,
    Poseidon,
    RescuePrime,
}

#[derive(Copy, Clone, Debug, serde::Serialize, serde::Deserialize)]
pub enum CustomGate {
    QuinticWidth4,
    QuinticWidth3,
    None,
}

#[derive(Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Sbox {
    Alpha(u64),
    AlphaInverse([u64; 4]), // TODO
}

impl std::fmt::Debug for Sbox {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Alpha(_) => write!(f, "quintic sbox"),
            Self::AlphaInverse(_) => write!(f, "quintic inverse sbox"),
        }
    }
}

pub trait HashParams<E: Engine, const RATE: usize, const WIDTH: usize>:
    Clone + Send + Sync + serde::Serialize + serde::de::DeserializeOwned
{
    fn hash_family(&self) -> HashFamily;
    fn constants_of_round(&self, round: usize) -> [E::Fr; WIDTH];
    fn mds_matrix(&self) -> [[E::Fr; WIDTH]; WIDTH];
    fn number_of_full_rounds(&self) -> usize;
    fn number_of_partial_rounds(&self) -> usize;
    fn alpha(&self) -> &Sbox;
    fn alpha_inv(&self) -> &Sbox;
    fn optimized_round_constants(&self) -> &[[E::Fr; WIDTH]];
    fn optimized_mds_matrixes(&self) -> (&[[E::Fr; WIDTH]; WIDTH], &[[[E::Fr; WIDTH]; WIDTH]]);
    fn custom_gate(&self) -> CustomGate;
    fn use_custom_gate(&mut self, gate: CustomGate);
}
