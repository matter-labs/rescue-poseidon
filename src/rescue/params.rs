use franklin_crypto::bellman::{Engine, PrimeField};

use crate::HasherParams;

pub fn rescue_params<E: Engine, const RATE: usize, const STATE_WIDTH: usize>() -> (HasherParams<E, RATE, STATE_WIDTH>, E::Fr, Option<E::Fr>) {
    let full_rounds = 22;
    let security_level = 126;

    let mut params = HasherParams::new(        
        security_level,
        full_rounds,
        0,
    );

    // let number_of_round_constants = (1 + 2 * full_rounds) * STATE_WIDTH;
    let rounds_tag = b"Rescue_f";
    let _mds_tag = b"ResM0003";
    params.compute_round_constants(full_rounds, rounds_tag);
    params.compute_mds_matrix_for_rescue();

    let alpha_u64 = 5u64;
    let alpha = E::Fr::from_str(&alpha_u64.to_string()).unwrap();
    let alpha_inv = crate::common::utils::compute_gcd::<E>(alpha_u64);

    (params, alpha, alpha_inv)
}


