use super::traits::SpongeGadget;
use crate::circuit::sponge::GenericSpongeGadget;
use crate::poseidon::PoseidonParams;
use crate::rescue::RescueParams;
use crate::rescue_prime::RescuePrimeParams;
use crate::sponge::GenericSponge;
use crate::tests::init_cs;
use crate::tests::init_rng;
use crate::traits::HashParams;
use crate::traits::Sponge;
use franklin_crypto::bellman::pairing::bn256::Bn256;
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem,
    bellman::{Engine, Field},
    plonk::circuit::allocated_num::AllocatedNum,
    plonk::circuit::allocated_num::Num,
};
use rand::{Rand, Rng};

fn test_sponge_gadget<
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
    let mut inputs = vec![E::Fr::zero(); 2];
    let mut inputs_as_num = vec![Num::Constant(E::Fr::zero()); 2];
    for (i1, i2) in inputs.iter_mut().zip(inputs_as_num.iter_mut()) {
        *i1 = E::Fr::rand(rng);
        *i2 = Num::Variable(AllocatedNum::alloc(cs, || Ok(*i1)).unwrap());
    }

    let mut hasher = GenericSponge::<_, P, RATE, WIDTH>::from_params(params);
    hasher.absorb(&inputs);
    let expected = hasher.squeeze(None);

    let mut gadget = GenericSpongeGadget::<_, P, RATE, WIDTH>::from(params);
    gadget.absorb(cs, &inputs_as_num).unwrap();
    let actual = gadget.squeeze(cs, None).unwrap();

    // cs.finalize();
    // assert!(cs.is_satisfied());
    // println!("last step number {}", cs.get_current_step_number());

    for (sponge, gadget) in expected.iter().zip(actual.iter()) {
        assert_eq!(gadget.get_value().unwrap(), *sponge);
    }
}

#[test]
fn test_sponge_gadget_rescue() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = RescueParams::default();
    test_sponge_gadget::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}

#[test]
fn test_sponge_gadget_poseidon() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = PoseidonParams::default();
    test_sponge_gadget::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}

#[test]
#[should_panic] // TODO
fn test_sponge_gadget_rescue_prime() {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let rng = &mut init_rng();
    let cs = &mut init_cs::<Bn256>();

    let params = RescuePrimeParams::default();
    test_sponge_gadget::<_, _, _, _, RATE, WIDTH>(cs, rng, &params);
}
