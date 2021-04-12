use franklin_crypto::bellman::bn256::{Bn256, Fr};
use franklin_crypto::bellman::Engine;
use rand::{Rand, SeedableRng, XorShiftRng};
#[allow(dead_code)]
use rescue_poseidon::generic_hash;
use rescue_poseidon::poseidon::PoseidonParams;
use rescue_poseidon::rescue::RescueParams;
use rescue_poseidon::rescue::{
    generic_rescue_hash, generic_rescue_var_length, rescue_hash, rescue_hash_var_length,
};
use rescue_poseidon::rescue_prime::RescuePrimeParams;
use rescue_poseidon::GenericSponge;
use rescue_poseidon::Sponge;
use std::convert::TryInto;

pub(crate) fn init_rng() -> XorShiftRng {
    const TEST_SEED: [u32; 4] = [0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654];
    XorShiftRng::from_seed(TEST_SEED)
}

fn main() {
    run_generic_hash_fixed_length::<Bn256>();
    run_rescue_fixed_length_example::<Bn256>();
    run_rescue_var_length_example::<Bn256>();
    run_generic_rescue_fixed_length_example::<Bn256>();
    run_generic_rescue_var_length_example::<Bn256>();
    run_generic_sponge_with_rescue_params::<Bn256>();
    run_generic_sponge_with_single_squeeze::<Bn256>();
    run_generic_sponge_with_requested_nuumber_output::<Bn256>();
}

fn run_generic_hash_fixed_length<E: Engine>() {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 2;
    let rng = &mut init_rng();
    let inputs: [Fr; INPUT_LENGTH] = (0..INPUT_LENGTH)
        .map(|_| Fr::rand(rng))
        .collect::<Vec<Fr>>()
        .try_into()
        .expect("constant array");
    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&rescue_params, &inputs);
    assert_eq!(result.len(), RATE);

    // now, hash with poseidon params
    let poseidon_params = PoseidonParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&poseidon_params, &inputs);
    assert_eq!(result.len(), RATE);

    // // now, hash with rescue prime params
    let rescue_prime_params = RescuePrimeParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&rescue_prime_params, &inputs);
    assert_eq!(result.len(), RATE);
}

fn run_rescue_fixed_length_example<E: Engine>() {
    const INPUT_LENGTH: usize = 2;
    let rng = &mut init_rng();
    let input = (0..INPUT_LENGTH)
        .map(|_| Fr::rand(rng))
        .collect::<Vec<Fr>>();

    let result = rescue_hash::<Bn256, INPUT_LENGTH>(&input.try_into().expect("static vector"));
    assert_eq!(result.len(), 2);
}

fn run_rescue_var_length_example<E: Engine>() {
    let rng = &mut init_rng();
    let input = (0..4).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();

    let result = rescue_hash_var_length::<Bn256>(&input);
    assert_eq!(result.len(), 2);
}

fn run_generic_rescue_fixed_length_example<E: Engine>() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    const INPUT_LENGTH: usize = 5;
    let rng = &mut init_rng();
    let input = (0..INPUT_LENGTH)
        .map(|_| Fr::rand(rng))
        .collect::<Vec<Fr>>();

    let result = generic_rescue_hash::<Bn256, RATE, WIDTH, INPUT_LENGTH>(
        &input.try_into().expect("static vector"),
    );
    assert_eq!(result.len(), RATE);
}

fn run_generic_rescue_var_length_example<E: Engine>() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    const INPUT_LENGTH: usize = 8; // input length should be multiple of RATE
    let rng = &mut init_rng();
    let input = (0..INPUT_LENGTH)
        .map(|_| Fr::rand(rng))
        .collect::<Vec<Fr>>();

    let result = generic_rescue_var_length::<Bn256, RATE, WIDTH>(&input);
    assert_eq!(result.len(), RATE);
}

fn run_generic_sponge_with_rescue_params<E: Engine>() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let rng = &mut init_rng();
    let input = (0..2).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();

    let new_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let mut hasher = GenericSponge::from_params(&new_params);
    hasher.absorb(&input);
    let result = hasher.squeeze(None);
    assert_eq!(result.len(), RATE);
}

fn run_generic_sponge_with_single_squeeze<E: Engine>() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let rng = &mut init_rng();
    let input = (0..2).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();

    let new_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let mut hasher = GenericSponge::from_params(&new_params);
    hasher.absorb(&input);
    let result = hasher.squeeze(Some(1));
    // Specifying output length may cause to lose some bits of hash result
    assert_eq!(result.len(), 1);
}

fn run_generic_sponge_with_requested_nuumber_output<E: Engine>() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let rng = &mut init_rng();
    let input = (0..2).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>();

    let requested_number_of_output: usize = 6;
    let new_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let mut hasher = GenericSponge::from_params(&new_params);
    hasher.absorb(&input);
    let result = hasher.squeeze(Some(requested_number_of_output));
    assert_eq!(result.len(), requested_number_of_output);
}
