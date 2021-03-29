use franklin_crypto::bellman::{Engine, Field, PrimeField};

use crate::HasherParams;

use crate::common::matrix::{compute_optimized_matrixes, mmul_assign, try_inverse};

pub fn poseidon_params<E: Engine, const STATE_WIDTH: usize, const RATE: usize>(
) -> (HasherParams<E,STATE_WIDTH, RATE>, E::Fr) {
    let security_level = 80;
    let full_rounds = 8;
    // let partial_rounds = 83;
    let partial_rounds = 33;

    let mut params = HasherParams::new(security_level, full_rounds, partial_rounds);

    let number_of_rounds = full_rounds + partial_rounds;
    let rounds_tag = b"Rescue_f";
    params.compute_round_constants(number_of_rounds, rounds_tag);
    params.compute_mds_matrix_for_poseidon();

    let alpha = E::Fr::from_str("5").unwrap();

    (params, alpha)
}

pub(crate) fn poseidon_light_params<E: Engine, const RATE: usize, const STATE_WIDTH: usize>() -> (
    HasherParams<E, STATE_WIDTH, RATE>,
    E::Fr,
    Vec<[E::Fr; STATE_WIDTH]>,
    (
        [[E::Fr; STATE_WIDTH]; STATE_WIDTH],
        Vec<[[E::Fr; STATE_WIDTH]; STATE_WIDTH]>,
    ),
) {
    let (params, alpha) = poseidon_params();

    let optimized_constants = compute_optimized_round_constants::<E, STATE_WIDTH>(
        params.round_constants(),
        &params.mds_matrix,
        params.partial_rounds,
        params.full_rounds,
    );

    // TODO:
    const SUBDIM: usize = 2;
    let optimized_matrixes = compute_optimized_matrixes::<E, STATE_WIDTH, SUBDIM>(
        params.partial_rounds,
        &params.mds_matrix,
    );
    (params, alpha, optimized_constants, optimized_matrixes)
}

// start from last round and walk to first round
// compute equivalent eq_k_i = MC^-1*k_i
// split it into two parts one for non-linear other for accumulation
// move it further to top
pub(crate) fn compute_optimized_round_constants<E: Engine, const STATE_WIDTH: usize>(
    constants: &[[E::Fr; STATE_WIDTH]],
    original_mds: &[[E::Fr; STATE_WIDTH]; STATE_WIDTH],
    number_of_partial_rounds: usize,
    number_of_full_rounds: usize,
) -> Vec<[E::Fr; STATE_WIDTH]> {
    let mds_inverse = try_inverse::<E, STATE_WIDTH>(original_mds).expect("has inverse");
    let number_of_half_rounds = number_of_full_rounds / 2;
    let start = number_of_half_rounds;
    let end = start + number_of_partial_rounds - 1;
    let mut acc: [E::Fr; STATE_WIDTH] = constants[end];
    let mut optimized_constants: Vec<[E::Fr; STATE_WIDTH]> = vec![];
    for round in (start..end).rev() {
        let mut inv = acc;
        mmul_assign::<E, STATE_WIDTH>(&mds_inverse, &mut inv);
        // make it two parts

        let mut second = [E::Fr::zero(); STATE_WIDTH];
        second[0] = inv[0];
        optimized_constants.push(second);

        let mut first = inv;
        first[0] = E::Fr::zero();

        // vector addition
        acc = [E::Fr::zero(); STATE_WIDTH];
        constants[round]
            .iter()
            .enumerate()
            .zip(first.iter())
            .for_each(|((idx, a), b)| {
                let mut tmp = a.clone();
                tmp.add_assign(&b);
                acc[idx] = tmp;
            });
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

    // let mut optimized_constants: Vec<[E::Fr; STATE_WIDTH]> =
    //     Vec::with_capacity(number_of_partial_rounds + number_of_full_rounds);

    // final_constants
    //     .chunks_exact(STATE_WIDTH)
    //     .zip(optimized_constants.iter_mut())
    //     .for_each(|(values, constants)| {
    //         *constants = values.try_into().expect("round constants in const");
    //     });

    final_constants
}
