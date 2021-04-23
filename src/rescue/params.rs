use franklin_crypto::bellman::{Engine, PrimeField};

use crate::common::params::InnerHashParameters;
use crate::traits::{HashParams, HashFamily};
use std::convert::TryInto;


#[derive(Clone, Debug)]
pub struct RescueParams<E: Engine, const RATE: usize, const WIDTH: usize> {
    pub(crate) full_rounds: usize,
    pub(crate) round_constants: Vec<[E::Fr; WIDTH]>,
    pub(crate) mds_matrix: [[E::Fr; WIDTH]; WIDTH],
    pub(crate) alpha: E::Fr,
    pub(crate) alpha_inv: E::Fr,
    pub(crate) allow_custom_gate: bool,
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
        let (params, alpha, alpha_inv) = compute_params::<E, RATE, WIDTH>();
        Self {
            full_rounds: params.full_rounds,
            round_constants: params
                .round_constants()
                .try_into()
                .expect("round constants"),
            mds_matrix: *params.mds_matrix(),
            alpha,
            alpha_inv,
            allow_custom_gate: true,
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

    fn alpha(&self) -> E::Fr {
        self.alpha
    }

    fn alpha_inv(&self) -> E::Fr {
        self.alpha_inv
    }

    fn optimized_mds_matrixes(&self) -> (&[[E::Fr; WIDTH]; WIDTH], &[[[E::Fr; WIDTH];WIDTH]]) {
        unimplemented!("Rescue doesn't use optimized matrixes")
    }

    fn optimized_round_constants(&self) -> &[[E::Fr; WIDTH]] {
        unimplemented!("Rescue doesn't use optimized round constants")
    }

    fn can_use_custom_gates(&self) -> bool {
        true
    }
    
    fn set_allow_custom_gate(&mut self, allow: bool) {
        self.allow_custom_gate = allow;    
    }
}


pub(crate) fn compute_params<E: Engine, const RATE: usize, const WIDTH: usize>() -> (InnerHashParameters<E, RATE, WIDTH>, E::Fr, E::Fr) {
    let full_rounds = 22;
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

    let alpha_u64 = 5u64;
    let alpha = E::Fr::from_str(&alpha_u64.to_string()).unwrap();
    let alpha_inv = crate::common::utils::compute_gcd::<E>(alpha_u64).expect("inverse of alpha");

    (params, alpha, alpha_inv)
}


