use crate::poseidon::params::PoseidonParams;
use crate::rescue::params::RescueParams;
use crate::rescue_prime::params::RescuePrimeParams;
use crate::tests::init_cs;
use crate::tests::init_rng;
use crate::traits::HashParams;
use franklin_crypto::bellman::pairing::bn256::Bn256;
use franklin_crypto::{bellman::plonk::better_better_cs::cs::ConstraintSystem, bellman::Engine};
use rand::{Rng, Rand};
use franklin_crypto::plonk::circuit::allocated_num::Num;
use franklin_crypto::bellman::Field;
use franklin_crypto::plonk::circuit::allocated_num::AllocatedNum;
use crate::circuit::sponge::CircuitGenericSponge;
use crate::sponge::GenericSponge;

fn test_circuit_var_len_generic_hasher<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, RATE, WIDTH>,
    R: Rng,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    rng: &mut R,
    params: &P,
) {
    let mut inputs = [E::Fr::zero(); 2];
    let mut inputs_as_num = [Num::Constant(E::Fr::zero()); 2];
    for (i1, i2) in inputs.iter_mut().zip(inputs_as_num.iter_mut()) {
        *i1 = E::Fr::rand(rng);
        *i2 = Num::Variable(AllocatedNum::alloc(cs, || Ok(*i1)).unwrap());
    }

    let mut hasher = GenericSponge::<_, P, RATE, WIDTH>::new_from_params(params);
    hasher.absorb_multiple(&inputs);
    let expected = hasher.squeeze().expect("a squeezed elem");

    let mut circuit_gadget = CircuitGenericSponge::<_, P, RATE, WIDTH>::new_from_params(params);
    circuit_gadget.absorb_multiple(cs, &inputs_as_num).unwrap();
    let actual = circuit_gadget.squeeze(cs).unwrap().expect("a squeezed elem");

    // cs.finalize();
    // assert!(cs.is_satisfied());
    // println!("last step number {}", cs.get_current_step_number());

    assert_eq!(actual.get_value().unwrap(), expected);

}
fn test_circuit_fixed_len_generic_hasher<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, RATE, WIDTH>,
    R: Rng,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    rng: &mut R,
    params: &P,
) {
    let mut inputs = [E::Fr::zero(); 2];
    let mut inputs_as_num = [Num::Constant(E::Fr::zero()); 2];
    for (i1, i2) in inputs.iter_mut().zip(inputs_as_num.iter_mut()) {
        *i1 = E::Fr::rand(rng);
        *i2 = Num::Variable(AllocatedNum::alloc(cs, || Ok(*i1)).unwrap());
    }

    let expected = GenericSponge::<_, P, RATE, WIDTH>::hash(&inputs, params);

    let actual = CircuitGenericSponge::<_, P, RATE, WIDTH>::hash(cs, &inputs_as_num, &params).unwrap();

    // cs.finalize();
    // assert!(cs.is_satisfied());
    // println!("last step number {}", cs.get_current_step_number());

    assert_eq!(actual[0].get_value().unwrap(), expected[0]);

}

#[test]
fn test_circuit_fixed_len_rescue_hasher() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = RescueParams::default();
    test_circuit_fixed_len_generic_hasher::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}

#[test]
fn test_circuit_fixed_len_poseidon_hasher() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = PoseidonParams::default();
    test_circuit_fixed_len_generic_hasher::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}


#[test]
fn test_circuit_fixed_len_rescue_prime_hasher() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = RescuePrimeParams::default();
    test_circuit_fixed_len_generic_hasher::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}
#[test]
fn test_circuit_var_len_rescue_hasher() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = RescueParams::default();
    test_circuit_var_len_generic_hasher::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}

#[test]
fn test_circuit_var_len_poseidon_hasher() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = PoseidonParams::default();
    test_circuit_var_len_generic_hasher::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}


#[test]
fn test_circuit_var_len_rescue_prime_hasher() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = RescuePrimeParams::default();
    test_circuit_var_len_generic_hasher::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}
