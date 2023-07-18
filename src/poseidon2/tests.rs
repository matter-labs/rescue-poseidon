use boojum::field::goldilocks::GoldilocksField;
use franklin_crypto::bellman::{pairing::bn256::{Bn256, Fr}, plonk::commitments};
use boojum::algebraic_props::round_function::AbsorptionModeTrait;
use boojum::field::SmallField;
use boojum::field::U64Representable;
use rand::Rand;
use rand::Rng;
use boojum::field::rand_from_rng;

use super::Poseidon2Sponge;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct TestingAbsorption;

impl AbsorptionModeTrait<Fr> for TestingAbsorption {
    #[inline(always)]
    fn absorb(dst: &mut Fr, src: &Fr) {
        *dst = *src;
    }
    #[inline(always)]
    fn pad(_dst: &mut Fr) {}
}

#[test]
fn test_different_absorbtions() {
    let num_elements = 1000;
    let mut rng = rand::thread_rng();
    let buffer: Vec<_> = (0..num_elements).map(|_| Fr::rand(&mut rng)).collect();

    // absorb by 1
    let mut hash = Poseidon2Sponge::<Bn256, GoldilocksField, TestingAbsorption, 2, 3>::new();
    for element in buffer.iter() {
        hash.absorb_single(element);
    }
    let commitment1 = hash.finalize();

    // absorb simultaneously
    let mut hash = Poseidon2Sponge::<Bn256, GoldilocksField, TestingAbsorption, 2, 3>::new();
    hash.absorb(&buffer);
    let commitment2 = hash.finalize();

    assert_eq!(commitment1, commitment2);
    dbg!(commitment1);
}

#[test]
fn test_vs_poseidon() {
    use crate::poseidon::{poseidon_hash, poseidon_round_function};
    use crate::poseidon2::{poseidon2_hash, poseidon2_round_function};

    const NUM_ELEMENTS: usize = 10000;
    let mut rng = rand::thread_rng();
    let buffer = [0; NUM_ELEMENTS].map(|_| Fr::rand(&mut rng));

    // hash by poseidon
    let start = std::time::Instant::now();
    poseidon_hash::<Bn256, NUM_ELEMENTS>(&buffer);
    dbg!(start.elapsed());

    // hash by poseidon2 
    // 61% better performance
    let start = std::time::Instant::now();
    poseidon2_hash::<Bn256, 3, 2, NUM_ELEMENTS>(&buffer);
    dbg!(start.elapsed());

    let mut buffer = buffer.to_vec();

    // round functions by poseidon
    let params = crate::PoseidonParams::<Bn256, 2, 3>::default();

    let start = std::time::Instant::now();
    for chunk in buffer.chunks_exact_mut(3) {
        poseidon_round_function::<Bn256, _, 2, 3>(&params, chunk.try_into().unwrap());
    }
    dbg!(start.elapsed());

    // round functions by poseidon2
    // 51% better performance
    let params = crate::poseidon2::Poseidon2Params::<Bn256, 2, 3>::default();

    let start = std::time::Instant::now();
    for chunk in buffer.chunks_exact_mut(3) {
        poseidon2_round_function::<Bn256, 2, 3>(chunk.try_into().unwrap(), &params);
    }
    dbg!(start.elapsed());
}

#[test]
fn test_of_sponge_state() {
    let num_elements = 5;
    let mut rng = rand::thread_rng();
    let buffer1: Vec<_> = (0..num_elements).map(|_| Fr::rand(&mut rng)).collect();

    let mut rng = rand::thread_rng();
    let buffer2: Vec<_> = (0..num_elements).map(|_| 
        GoldilocksField::from_u64_unchecked(rng.gen_range(0, GoldilocksField::CHAR))
    ).collect();

    dbg!(&buffer1, &buffer2);

    let mut hash = Poseidon2Sponge::<Bn256, GoldilocksField, TestingAbsorption, 2, 3>::new();

    dbg!(&hash);

    hash.absorb_single(&buffer1[0]);
    dbg!(&hash);

    hash.absorb_single_small_field(&buffer2[0]);
    dbg!(&hash);

    hash.absorb_single_small_field(&buffer2[1]);
    dbg!(&hash);

    hash.absorb_single(&buffer1[1]);
    dbg!(&hash);

    hash.absorb_single(&buffer1[2]);
    dbg!(&hash);

    hash.absorb_single_small_field(&buffer2[2]);
    dbg!(&hash);

    dbg!(&hash.finalize());
    dbg!(&hash);
}
