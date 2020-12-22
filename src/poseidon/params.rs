use franklin_crypto::bellman::{Engine, Field, PrimeField};

use crate::HasherParams;

use crate::common::matrix::{compute_optimized_matrixes, mmul_assign, try_inverse};

pub fn poseidon_params<E: Engine>() -> (HasherParams<E>, E::Fr) {
    let rate = 2;
    let capacity = 1;
    let state_width = rate + capacity;
    let security_level = 80;
    let full_rounds = 8;
    // let partial_rounds = 83;
    let partial_rounds = 33;

    let mut params = HasherParams::new(rate, capacity, security_level, full_rounds, partial_rounds);

    let number_of_round_constants = (full_rounds + partial_rounds) * state_width;
    let rounds_tag = b"Rescue_f";
    params.compute_round_constants(number_of_round_constants, rounds_tag);
    params.compute_mds_matrix_for_poseidon();

    let alpha = E::Fr::from_str("5").unwrap();

    (params, alpha)
}

pub(crate) fn poseidon_light_params<E: Engine>() -> (
    HasherParams<E>,
    E::Fr,
    Vec<Vec<E::Fr>>,
    (Vec<Vec<E::Fr>>, Vec<Vec<Vec<E::Fr>>>),
) {
    let (params, alpha) = poseidon_params();

    let optimized_constants = compute_optimized_round_constants::<E>(
        params.round_constants(),
        &params.mds_matrix,
        params.partial_rounds,
        params.full_rounds,
        params.state_width,
    );
    let optimized_matrixes =
        compute_optimized_matrixes::<E>(params.partial_rounds, &params.mds_matrix);
    (params, alpha, optimized_constants, optimized_matrixes)
}

// start from last round and walk to first round
// compute equivalent eq_k_i = MC^-1*k_i
// split it into two parts one for non-linear other for accumulation
// move it further to top
pub(crate) fn compute_optimized_round_constants<E: Engine>(
    constants: &[Vec<E::Fr>],
    original_mds: &[Vec<E::Fr>],
    number_of_partial_rounds: usize,
    number_of_full_rounds: usize,
    state_width: usize,
) -> Vec<Vec<E::Fr>> {
    let mds_inverse = try_inverse::<E>(original_mds).expect("has inverse");
    let number_of_half_rounds = number_of_full_rounds / 2;
    let start = number_of_half_rounds;
    let end = start + number_of_partial_rounds - 1;
    let mut acc = constants[end].to_vec();
    let mut optimized_constants = vec![];
    for round in (start..end).rev() {
        let mut inv = acc.to_vec();
        mmul_assign::<E>(&mds_inverse, &mut inv);
        // make it two parts

        let mut second = vec![E::Fr::zero(); state_width];
        second[0] = inv[0];
        optimized_constants.push(second);

        let mut first = inv;
        first[0] = E::Fr::zero();

        // vector addition
        acc = constants[round]
            .iter()
            .zip(first)
            .map(|(a, b)| {
                let mut tmp = a.clone();
                tmp.add_assign(&b);
                tmp
            })
            .collect();
    }
    optimized_constants.push(acc);
    optimized_constants.reverse();

    let mut final_constants = constants.to_vec();
    final_constants[start..end + 1]
        .iter_mut()
        .zip(optimized_constants)
        .for_each(|(a, b)| {
            *a = b;
        });

    final_constants
}
