use franklin_crypto::bellman::{Engine};

use crate::common::params::InnerHashParameters;
use crate::traits::{HashParams, HashFamily, Sbox, CustomGate};
use std::convert::TryInto;


#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct RescueParams<E: Engine, const RATE: usize, const WIDTH: usize> {
    pub(crate) full_rounds: usize,
    #[serde(serialize_with = "crate::serialize_vec_of_arrays")]
    #[serde(deserialize_with = "crate::deserialize_vec_of_arrays")]
    pub(crate) round_constants: Vec<[E::Fr; WIDTH]>,
    #[serde(serialize_with = "crate::serialize_array_of_arrays")]
    #[serde(deserialize_with = "crate::deserialize_array_of_arrays")]
    pub(crate) mds_matrix: [[E::Fr; WIDTH]; WIDTH],
    pub(crate) alpha: Sbox,
    pub(crate) alpha_inv: Sbox,
    pub(crate) custom_gate: CustomGate,
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> PartialEq for RescueParams<E, RATE, WIDTH>{
    fn eq(&self, other: &Self) -> bool {
        self.hash_family() == other.hash_family()
    }
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> Default
    for RescueParams<E, RATE, WIDTH>
{
    fn default() -> Self {
        let (params, alpha, alpha_inv) = compute_params::<E, RATE, WIDTH, 4>();
        Self {
            full_rounds: params.full_rounds,
            round_constants: params
                .round_constants()
                .try_into()
                .expect("round constants"),
            mds_matrix: *params.mds_matrix(),
            alpha: Sbox::Alpha(alpha),
            alpha_inv: Sbox::AlphaInverse(alpha_inv),
            custom_gate: CustomGate::None,
        }
    }
}


impl<E: Engine, const RATE: usize, const WIDTH: usize> HashParams<E, RATE, WIDTH>
    for RescueParams<E, RATE, WIDTH>
{
    fn hash_family(&self) -> HashFamily {
        HashFamily::Rescue
    }

    fn constants_of_round(&self, round: usize) -> [E::Fr; WIDTH] {
        self.round_constants[round]
    }

    fn mds_matrix(&self) -> [[E::Fr; WIDTH]; WIDTH] {
        self.mds_matrix
    }

    fn number_of_full_rounds(&self) -> usize {
        self.full_rounds
    }

    fn number_of_partial_rounds(&self) -> usize {
        unimplemented!("Rescue doesn't have partial rounds.")
    }

    fn alpha(&self) -> &Sbox {
        &self.alpha
    }

    fn alpha_inv(&self) -> &Sbox {
        &self.alpha_inv
    }

    fn optimized_mds_matrixes(&self) -> (&[[E::Fr; WIDTH]; WIDTH], &[[[E::Fr; WIDTH];WIDTH]]) {
        unimplemented!("Rescue doesn't use optimized matrixes")
    }

    fn optimized_round_constants(&self) -> &[[E::Fr; WIDTH]] {
        unimplemented!("Rescue doesn't use optimized round constants")
    }

    fn custom_gate(&self) -> CustomGate {
        self.custom_gate
    }
    
    fn use_custom_gate(&mut self, custom_gate: CustomGate) {
        self.custom_gate = custom_gate;    
    }
}


pub(crate) fn compute_params<E: Engine, const RATE: usize, const WIDTH: usize, const N: usize>() -> (InnerHashParameters<E, RATE, WIDTH>, u64, [u64; N]) {
    // let full_rounds = 22;
    let full_rounds = 8;
    let security_level = 126;

    let mut params = InnerHashParameters::new(        
        security_level,
        full_rounds,
        0,
    );

    let rounds_tag = b"Rescue_f";
    let _mds_tag = b"ResM0003";
    let total_number_of_rounds = 2*full_rounds + 1;
    
    params.compute_round_constants(total_number_of_rounds, rounds_tag);
    params.compute_mds_matrix_for_rescue();

    let alpha = 5u64;
    let alpha_inv = crate::common::utils::compute_gcd::<E, N>(alpha).expect("inverse of alpha");

    (params, alpha, alpha_inv)
}


