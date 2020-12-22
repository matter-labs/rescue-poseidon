use franklin_crypto::bellman::{Engine, PrimeField};

use crate::HasherParams;

pub fn rescue_params<E: Engine>() -> (HasherParams<E>, E::Fr, Option<E::Fr>) {
    let rate = 2;
    let capacity = 1;
    let state_width = rate + capacity;
    let full_rounds = 22;
    let security_level = 126;

    let mut params = HasherParams::new(
        rate,
        capacity,
        security_level,
        full_rounds,
        0,
    );

    let number_of_round_constants = (1 + 2 * full_rounds) * state_width;
    let rounds_tag = b"Rescue_f";
    let _mds_tag = b"ResM0003";
    params.compute_round_constants(number_of_round_constants, rounds_tag);
    params.compute_mds_matrix_for_rescue();

    let alpha_u64 = 5u64;
    let alpha = E::Fr::from_str(&alpha_u64.to_string()).unwrap();
    let alpha_inv = crate::common::utils::compute_gcd::<E>(alpha_u64);

    (params, alpha, alpha_inv)
}
