use franklin_crypto::bellman::{Engine, PrimeField};

use crate::common::params::HasherParams;

pub fn rescue_params<E: Engine, const RATE: usize, const WIDTH: usize>() -> (HasherParams<E, RATE, WIDTH>, E::Fr, E::Fr) {
    let full_rounds = 22;
    let security_level = 126;

    let mut params = HasherParams::new(        
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


