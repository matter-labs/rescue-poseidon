use crate::common::params::HasherParams;
use franklin_crypto::bellman::pairing::ff::{PrimeFieldRepr, ScalarEngine};
use franklin_crypto::bellman::pairing::Engine;
extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;
use franklin_crypto::bellman::pairing::bn256::Bn256;
use franklin_crypto::bellman::{Field, PrimeField};
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::{ExtendedGcd, Integer};
use num_traits::{One, Zero};
use std::ops::{Mul, Sub};
use std::convert::TryInto;

fn get_number_of_rounds(m: usize, r: usize, security_level: usize, alpha: usize) -> usize {
    let capacity = m - r;
    fn factorial(n: &BigUint) -> BigUint {
        if n.is_zero() {
            return BigUint::one();
        } else {
            return n * factorial(&(n - &BigUint::one()));
        }
    }

    fn binomial(n: &BigUint, k: &BigUint) -> BigUint {
        factorial(n) / (factorial(&(n - k)) * factorial(k))
    }

    let rate = m - capacity;
    let dcon = |n: usize| -> BigUint {
        let tmp = ((alpha - 1) * m * (n - 1)) as f64;
        let result = BigUint::from((0.5 * tmp).floor() as u64) + BigUint::from(2u8);
        BigUint::from(result)
    };

    let v = |n| -> BigUint {
        let result = m * (n - 1) + rate;
        BigUint::from(result)
    };

    let target = BigUint::from(2u128.pow(security_level as u32));

    let mut actual_l1 = 0;
    for l1 in 1..25 {
        if (binomial(&(v(l1) + dcon(l1)), &v(l1)).pow(2u32)) > target {
            actual_l1 = l1;
            break;
        }
    }
    assert!(actual_l1 > 0, "l1 must be greater than zero");

    (1.5 * actual_l1.max(5) as f64).ceil() as usize
}

fn compute_alpha(p: &[u8]) -> (BigUint, BigUint) {
    let p_big = BigInt::from_bytes_le(Sign::Plus, p);
    let p_minus_one = p_big.sub(BigInt::from(1));
    let mut actual_alpha = BigInt::from(0);
    for alpha in num_iter::range_inclusive(BigInt::from(3), p_minus_one.to_owned()) {
        if p_minus_one.gcd(&alpha).is_one() {
            actual_alpha = alpha;
            break;
        }
    }
    let ExtendedGcd {
        gcd,
        y: mut alpha_inv,
        ..
    } = p_minus_one.extended_gcd(&actual_alpha);
    assert!(gcd.is_one());
    if alpha_inv < BigInt::zero() {
        alpha_inv += p_minus_one;
    };

    (
        actual_alpha.to_biguint().unwrap(),
        alpha_inv.to_biguint().unwrap(),
    )
}

fn compute_round_constants<E: Engine, const STATE_WIDTH: usize>(
    modulus_bytes: &[u8],
    p_big: BigInt,
    m: usize,
    capacity: usize,
    security_level: usize,
    n: usize,
) -> Vec<[E::Fr; STATE_WIDTH]> {
    fn shake256(input: &[u8], num_bytes: usize) -> Box<[u8]> {
        use sha3::digest::ExtendableOutput;
        use sha3::digest::Update;
        use sha3::Shake256;

        let mut shake = Shake256::default();
        shake.update(input);
        shake.finalize_boxed(num_bytes)
    }

    let modulus_bit_len = (modulus_bytes.len() * 8 - 2) as f32;

    let bytes_per_int = ((modulus_bit_len / 8f32) + 1f32).ceil() as usize;
    let num_bytes = bytes_per_int * 2 * m * n;
    let seed_string = format!(
        "Rescue-XLIX({},{},{},{})",
        p_big, m, capacity, security_level
    );
    let seed_bytes = seed_string.as_bytes();
    let byte_string = shake256(seed_bytes, num_bytes);
    let mut round_constants = vec![];
    for i in 0..2 * m * n {
        let chunk = byte_string[bytes_per_int * i..bytes_per_int * (i + 1)].to_vec();
        let constant = chunk
            .iter()
            .enumerate()
            .map(|(i, j)| {
                // (256^j) * ZZ(chunk[j])
                let pow = BigUint::from(256u16).pow(i as u32);
                let zz = BigUint::from_bytes_le(&[*j]);
                let result = pow.clone().mul(zz.clone());

                result
            })
            .fold(BigUint::zero(), |acc, next| acc + next);
        let remainder = constant.mod_floor(&p_big.to_biguint().expect("valid modulus"));
        let mut bytes_le = remainder.to_bytes_le();
        if bytes_le.len() < 64 {
            bytes_le.resize(64, 0u8);
        }

        let mut repr = <E::Fr as PrimeField>::Repr::default();
        repr.read_le(&bytes_le[..]).unwrap();
        let constant_fe = E::Fr::from_repr(repr).unwrap();
        round_constants.push(constant_fe);
    }
    let mut final_constants = vec![[E::Fr::zero(); STATE_WIDTH]; n];

    round_constants
        .chunks_exact(STATE_WIDTH)
        .zip(final_constants.iter_mut())
        .for_each(|(src, dst)| *dst = src.try_into().expect("constants in const"));

    final_constants
}

pub fn rescue_prime_params<E: Engine, const STATE_WIDTH: usize, const RATE: usize>(
) -> (HasherParams<E, STATE_WIDTH, RATE>, E::Fr, E::Fr) {
    let security_level = 80;

    let mut modulus_bytes = vec![];
    let p_fe = <Bn256 as ScalarEngine>::Fr::char();
    p_fe.write_le(&mut modulus_bytes).unwrap();
    let p_big = BigInt::from_bytes_le(Sign::Plus, &modulus_bytes);
    let (alpha, alpha_inv) = compute_alpha(&modulus_bytes);
    let alpha = alpha.to_u32_digits()[0] as usize;
    let number_of_rounds = get_number_of_rounds(STATE_WIDTH, RATE, security_level, alpha);


    let mut params = HasherParams::new(security_level, number_of_rounds, 0);
    params.round_constants = compute_round_constants::<E, STATE_WIDTH>(
        &modulus_bytes,
        p_big,        
        STATE_WIDTH,
        STATE_WIDTH - RATE, // TODO: remove them
        security_level,
        number_of_rounds,
    );

    params.compute_mds_matrix_for_rescue();

    let alpha_fe = E::Fr::from_str(&alpha.to_string()).unwrap();
    let mut repr = <E::Fr as PrimeField>::Repr::default();
    repr.read_le(&alpha_inv.to_bytes_le()[..]).unwrap();
    let alpha_inv_fe = E::Fr::from_repr(repr).unwrap();

    (params, alpha_fe, alpha_inv_fe)
}

#[test]
fn test_calculate_number_of_rounds() {
    let p_fe = <Bn256 as ScalarEngine>::Fr::char();
    let m = 3;
    let capacity = 1;
    let security_level = 80;
    let mut modulus_bytes = vec![];
    p_fe.write_le(&mut modulus_bytes).unwrap();
    let p_big = BigInt::from_bytes_le(Sign::Plus, &modulus_bytes);
    let (alpha, alpha_inv) = compute_alpha(&modulus_bytes);
    let alpha = alpha.to_u32_digits()[0] as usize;
    let n = get_number_of_rounds(m, capacity, security_level, alpha);

    println!(
        "alpha {} alpha inv {:x} number of rounds {}",
        alpha, alpha_inv, n
    );
    let round_constants =
        compute_round_constants::<Bn256, 3>(&modulus_bytes, p_big, m, capacity, security_level, n);

    println!("number of rounds {}", n);
    println!("number of round constants {}", round_constants.len());
}
