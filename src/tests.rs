use crate::rescue::RescueHasher;
use crate::sponge::StatefulSponge;
use crate::{poseidon::PoseidonHasher, sponge::SpongeModes};
use franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
use franklin_crypto::bellman::{Field, PrimeField};
use franklin_crypto::group_hash::BlakeHasher;
use franklin_crypto::rescue::{bn256::Bn256RescueParams, RescueHashParams, StatefulRescue};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::{TrivialAssembly, Width4MainGateWithDNext},
    bellman::Engine,
    plonk::circuit::Width4WithCustomGates,
};
use poseidon_hash::StatefulSponge as PoseidonSponge;
use poseidon_hash::{
    bn256::Bn256PoseidonParams, PoseidonHashParams, StatefulSponge as StatefulSpongeOld,
};
use rand::{SeedableRng, XorShiftRng, Rand};

pub(crate) fn init_rng() -> XorShiftRng {
    XorShiftRng::from_seed(crate::common::TEST_SEED)
}
pub(crate) fn init_cs<E: Engine>(
) -> TrivialAssembly<E, Width4WithCustomGates, Width4MainGateWithDNext> {
    TrivialAssembly::<E, Width4WithCustomGates, Width4MainGateWithDNext>::new()
}
#[test]
fn test_rescue_bn256_fixed_length() {
    let rng = &mut init_rng();
    let input = (0..2).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();

    let old_params = Bn256RescueParams::new_checked_2_into_1();
    let expected = franklin_crypto::rescue::rescue_hash::<Bn256>(&old_params, &input);

    let actual = crate::rescue::rescue_hash::<Bn256>(&input);
    assert_eq!(expected[0], actual[0]);
}
#[test]
fn test_poseidon_bn256_fixed_length() {
    let rng = &mut init_rng();
    let input = (0..2).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();

    let old_params = Bn256PoseidonParams::new_checked_2_into_1();
    let expected = poseidon_hash::poseidon_hash::<Bn256>(&old_params, &input);

    let actual = crate::poseidon::poseidon_hash::<Bn256>(&input);
    assert_eq!(expected[0], actual[0]);
}

#[test]
fn test_rescue_params() {
    let old_params = Bn256RescueParams::new_checked_2_into_1();
    let (new_params, _, _) = crate::rescue::params::rescue_params::<Bn256>();

    let number_of_rounds = new_params.full_rounds;

    for round in 0..number_of_rounds {
        assert_eq!(
            old_params.round_constants(round as u32),
            new_params.constants_of_round(round)
        )
    }

    for row in 0..new_params.state_width {
        assert_eq!(
            old_params.mds_matrix_row(row as u32),
            new_params.mds_matrix[row]
        );
    }
}

#[test]
fn test_poseidon_params() {
    let old_params = Bn256PoseidonParams::new_checked_2_into_1();
    let (new_params, _) = crate::poseidon::params::poseidon_params::<Bn256>();

    let number_of_rounds = new_params.full_rounds;

    for round in 0..number_of_rounds {
        assert_eq!(
            old_params.round_constants(round as u32),
            new_params.constants_of_round(round)
        )
    }

    for row in 0..new_params.state_width {
        assert_eq!(
            old_params.mds_matrix_row(row as u32),
            new_params.mds_matrix[row]
        );
    }
}

#[test]
fn test_poseidon_comparisons_with_original_one() {
    let mut el = Fr::one();
    el.double();

    let input = vec![el; 2];

    let original_params = Bn256PoseidonParams::new_checked_2_into_1();
    let mut original_poseidon = PoseidonSponge::<Bn256>::new(&original_params);
    original_poseidon.absorb(&input);
    let expected = original_poseidon.squeeze_out_single();

    let mut hasher = PoseidonHasher::<Bn256>::default();
    hasher.absorb_multi(&input);
    let actual = hasher.squeeze();

    assert_eq!(actual[0], expected);
}

#[test]
fn test_rescue_comparisons_with_original_one() {
    let mut el = Fr::one();
    el.double();

    let input = vec![el; 2];

    let original_params = Bn256RescueParams::new_checked_2_into_1();
    let mut original_rescue = StatefulRescue::<Bn256>::new(&original_params);
    original_rescue.absorb(&input);
    let expected = original_rescue.squeeze_out_single();

    let mut hasher = RescueHasher::<Bn256>::default();
    hasher.absorb_multi(&input);
    let actual = hasher.squeeze();

    assert_eq!(actual[0], expected);
}

#[test]
fn test_poseidon_duplex() {
    use franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
    let mut el = Fr::one();
    el.double();

    let input = vec![el; 2];

    let original_params = Bn256PoseidonParams::new_checked_2_into_1();
    let mut original_poseidon = PoseidonSponge::<Bn256>::new(&original_params);
    original_poseidon.absorb_single_value(input[0]);
    let expected = original_poseidon.squeeze_out_single();

    let mut hasher = PoseidonHasher::<Bn256>::new(SpongeModes::Duplex(vec![]));
    hasher.absorb(input[0]);
    let actual = hasher.squeeze_single();

    assert_eq!(actual, expected);
}

#[test]
fn test_rescue_duplex() {
    use franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
    let mut el = Fr::one();
    el.double();

    let input = vec![el; 2];

    let original_params = Bn256RescueParams::new_checked_2_into_1();
    let mut original_rescue = StatefulRescue::<Bn256>::new(&original_params);
    original_rescue.absorb_single_value(input[0]);
    original_rescue.pad_if_necessary();
    let mut expected = vec![Fr::zero(); 2];
    expected[0] = original_rescue.squeeze_out_single();
    expected[1] = original_rescue.squeeze_out_single();

    let mut hasher = RescueHasher::<Bn256>::new(SpongeModes::Duplex(Vec::with_capacity(2)));
    hasher.absorb(input[0]);
    let mut actual = vec![Fr::zero(); 2];
    actual[0] = hasher.squeeze_single();
    actual[1] = hasher.squeeze_single();

    assert_eq!(actual, expected);
}
