use franklin_crypto::bellman::bn256::{Bn256, Fr};
use franklin_crypto::bellman::Engine;
use franklin_crypto::bellman::Field;
use franklin_crypto::plonk::circuit::allocated_num::{AllocatedNum, Num};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::{TrivialAssembly, Width4MainGateWithDNext},
    plonk::circuit::Width4WithCustomGates,
};
use rand::{Rand, SeedableRng, XorShiftRng};
use rescue_poseidon::PoseidonParams;
use rescue_poseidon::RescueParams;
use rescue_poseidon::RescuePrimeParams;
#[allow(dead_code)]
use rescue_poseidon::{generic_hash, CircuitGenericSponge, GenericSponge};
use std::convert::TryInto;

pub(crate) fn init_rng() -> XorShiftRng {
    const TEST_SEED: [u32; 4] = [0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654];
    XorShiftRng::from_seed(TEST_SEED)
}

pub(crate) fn init_cs<E: Engine>(
) -> TrivialAssembly<E, Width4WithCustomGates, Width4MainGateWithDNext> {
    TrivialAssembly::<E, Width4WithCustomGates, Width4MainGateWithDNext>::new()
}

fn test_input<E: Engine, const L: usize>() -> ([E::Fr; L], [Num<E>; L]) {
    let rng = &mut init_rng();
    let mut input = [E::Fr::zero(); L];
    let mut input_as_nums = [Num::Constant(E::Fr::zero()); L];
    for (inp, as_num) in input.iter_mut().zip(input_as_nums.iter_mut()) {
        *inp = E::Fr::rand(rng);
        *as_num = Num::Constant(*inp);
    }

    (input, input_as_nums)
}

fn main() {
    run_generic_hash_fixed_length::<Bn256>();
    run_generic_hash_var_length::<Bn256>();
    run_circuit_generic_hash_fixed_length::<Bn256>();
    run_circuit_generic_hash_var_length::<Bn256>();
}

fn run_simple_fixed_len_rescue_hash() {
    use franklin_crypto::bellman::bn256::Fr;
    use franklin_crypto::bellman::Field;
    use rescue_poseidon::rescue_hash;

    const INPUT_LENGTH: usize = 2;
    let input = [Fr::one(); INPUT_LENGTH]; // dummy input

    let result = rescue_hash::<Bn256, INPUT_LENGTH>(&input);
    assert_eq!(result.len(), 2);
}

fn run_generic_hash_fixed_length<E: Engine>() {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 2;
    let rng = &mut init_rng();
    let (input, _) = test_input::<Bn256, INPUT_LENGTH>();
    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&rescue_params, &input);
    assert_eq!(result.len(), RATE);

    // now, hash with poseidon params
    let poseidon_params = PoseidonParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&poseidon_params, &input);
    assert_eq!(result.len(), RATE);

    // now, hash with rescue prime params
    let rescue_prime_params = RescuePrimeParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&rescue_prime_params, &input);
    assert_eq!(result.len(), RATE);
}

fn run_generic_hash_var_length<E: Engine>() {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 4;
    let rng = &mut init_rng();
    let (input, _) = test_input::<Bn256, INPUT_LENGTH>();

    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let mut rescue_hasher = GenericSponge::new_from_params(&rescue_params);
    rescue_hasher.absorb_multiple(&input);
    let rescue_result = rescue_hasher.squeeze().expect("squeezed eleme");

    // go with poseidon
    let poseidon_params = PoseidonParams::<Bn256, RATE, WIDTH>::default();
    let mut poseidon_hasher = GenericSponge::new_from_params(&poseidon_params);
    poseidon_hasher.absorb_multiple(&input);
    let poseidon_result = poseidon_hasher.squeeze().expect("squeezed eleme");
}

fn run_circuit_generic_hash_fixed_length<E: Engine>() {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 2;
    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let (input, input_as_nums) = test_input::<E, INPUT_LENGTH>();
    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<E, RATE, WIDTH>::default();
    let result = GenericSponge::hash(&input, &rescue_params);
    assert_eq!(result.len(), RATE);

    // now, hash with poseidon params
    let poseidon_params = PoseidonParams::<E, RATE, WIDTH>::default();
    let result = GenericSponge::hash(&input, &poseidon_params);
    assert_eq!(result.len(), RATE);

    // now, hash with rescue prime params
    let rescue_prime_params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    let result = GenericSponge::hash(&input, &rescue_prime_params);
    assert_eq!(result.len(), RATE);
}

fn run_circuit_generic_hash_var_length<E: Engine>() {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 2;
    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let (input, input_as_nums) = test_input::<E, INPUT_LENGTH>();
    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<E, RATE, WIDTH>::default();
    let mut rescue_hasher = GenericSponge::new_from_params(&rescue_params);
    for inp in input.iter() {
        rescue_hasher.absorb(*inp);
    }
    let _ = rescue_hasher.squeeze().unwrap();

    // now, hash with poseidon params
    let poseidon_params = PoseidonParams::<E, RATE, WIDTH>::default();
    let mut poseidon_hasher = GenericSponge::new_from_params(&poseidon_params);
    for inp in input.iter() {
        poseidon_hasher.absorb(*inp);
    }
    let _ = poseidon_hasher.squeeze().unwrap();

    // now, hash with rescue prime params
    let rescue_prime_params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    let mut rescue_prime_hasher = GenericSponge::new_from_params(&rescue_prime_params);
    for inp in input.iter() {
        rescue_prime_hasher.absorb(*inp);
    }
    let _ = rescue_prime_hasher.squeeze().unwrap();
}
