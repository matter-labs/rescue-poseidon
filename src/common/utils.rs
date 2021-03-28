use franklin_crypto::bellman::pairing::ff::{Field, PrimeField};
use franklin_crypto::bellman::Engine;
use rand::Rng;
extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;

use self::num_bigint::{BigInt, BigUint};
use self::num_integer::{ExtendedGcd, Integer};
use self::num_traits::{One, ToPrimitive, Zero};
use std::convert::TryInto;

// Batch inverses vector of elements required for MDS matrix.
pub(crate) fn batch_inversion<E: Engine>(v: &mut [E::Fr]) {
    // Montgomery’s Trick and Fast Implementation of Masked AES
    // Genelle, Prouff and Quisquater
    // Section 3.2

    // First pass: compute [a, ab, abc, ...]
    let mut prod = Vec::with_capacity(v.len());
    let mut tmp = E::Fr::one();
    for g in v
        .iter()
        // Ignore zero elements
        .filter(|g| !g.is_zero())
    {
        tmp.mul_assign(&g);
        prod.push(tmp);
    }

    // Invert `tmp`.
    tmp = tmp.inverse().unwrap(); // Guaranteed to be nonzero.

    // Second pass: iterate backwards to compute inverses
    for (g, s) in v
        .iter_mut()
        // Backwards
        .rev()
        // Ignore normalized elements
        .filter(|g| !g.is_zero())
        // Backwards, skip last element, fill in one for last term.
        .zip(prod.into_iter().rev().skip(1).chain(Some(E::Fr::one())))
    {
        // tmp := tmp * g.z; g.z := tmp * s = 1/z
        let mut newtmp = tmp;
        newtmp.mul_assign(&g);
        *g = tmp;
        g.mul_assign(&s);
        tmp = newtmp;
    }
}

// Computes scalar product of two same length vector.
pub(crate) fn scalar_product<E: Engine>(a: &[E::Fr], b: &[E::Fr]) -> E::Fr {
    let mut acc = E::Fr::zero();
    for (a, b) in a.iter().zip(b.iter()) {
        let mut tmp = a.clone();
        tmp.mul_assign(&b);
        acc.add_assign(&tmp);
    }
    acc
}

// Construct MDS matrix which required by lineary layer of permutation function.
pub(crate) fn construct_mds_matrix<E: Engine, R: Rng, const S: usize>(
    rng: &mut R,
) -> [[E::Fr; S]; S] {
    let state_width = S;

    loop {
        let x: Vec<E::Fr> = (0..state_width).map(|_| rng.gen()).collect();
        let y: Vec<E::Fr> = (0..state_width).map(|_| rng.gen()).collect();

        let mut invalid = false;

        // quick and dirty check for uniqueness of x
        for i in 0..(state_width) {
            if invalid {
                continue;
            }
            let el = x[i];
            for other in x[(i + 1)..].iter() {
                if el == *other {
                    invalid = true;
                    break;
                }
            }
        }

        if invalid {
            continue;
        }

        // quick and dirty check for uniqueness of y
        for i in 0..(state_width) {
            if invalid {
                continue;
            }
            let el = y[i];
            for other in y[(i + 1)..].iter() {
                if el == *other {
                    invalid = true;
                    break;
                }
            }
        }

        if invalid {
            continue;
        }

        // quick and dirty check for uniqueness of x vs y
        for i in 0..(state_width) {
            if invalid {
                continue;
            }
            let el = x[i];
            for other in y.iter() {
                if el == *other {
                    invalid = true;
                    break;
                }
            }
        }

        if invalid {
            continue;
        }

        // by previous checks we can be sure in uniqueness and perform subtractions easily
        let mut mds_matrix = vec![E::Fr::zero(); state_width * state_width];
        for (i, x) in x.into_iter().enumerate() {
            for (j, y) in y.iter().enumerate() {
                let place_into = i * (state_width) + j;
                let mut element = x;
                element.sub_assign(y);
                mds_matrix[place_into] = element;
            }
        }

        // now we need to do the inverse
        batch_inversion::<E>(&mut mds_matrix[..]);

        let mut result = [[E::Fr::zero(); S]; S];

        mds_matrix
            .chunks_exact(S)
            .zip(result.iter_mut())
            .for_each(|(values, row)| *row = values.try_into().expect("row in const"));

        return result;
    }
}

// Computes GCD of an element. It basically computes inverse of alpha in given finite field.
pub(crate) fn compute_gcd<E: Engine>(n: u64) -> Option<E::Fr> {
    let n_big = BigUint::from(n);

    let mut p_minus_one_biguint = BigUint::from(0u64);
    for limb in E::Fr::char().as_ref().iter().rev() {
        p_minus_one_biguint <<= 64;
        p_minus_one_biguint += BigUint::from(*limb);
    }

    p_minus_one_biguint -= BigUint::one();

    fn biguint_to_u64_array(mut v: BigUint) -> [u64; 4] {
        let m: BigUint = BigUint::from(1u64) << 64;
        let mut ret = [0; 4];

        for idx in 0..4 {
            ret[idx] = (&v % &m).to_u64().unwrap();
            v >>= 64;
        }
        assert!(v.is_zero());
        ret
    }

    let alpha_signed = BigInt::from(n_big);
    let p_minus_one_signed = BigInt::from(p_minus_one_biguint);

    let ExtendedGcd { gcd, x: _, y, .. } = p_minus_one_signed.extended_gcd(&alpha_signed);
    assert!(gcd.is_one());
    let y = if y < BigInt::zero() {
        let mut y = y;
        y += p_minus_one_signed;

        y.to_biguint().expect("must be > 0")
    } else {
        y.to_biguint().expect("must be > 0")
    };

    let inv_alpha = biguint_to_u64_array(y);

    let mut alpha_inv_repr = <E::Fr as PrimeField>::Repr::default();
    for (r, limb) in alpha_inv_repr.as_mut().iter_mut().zip(inv_alpha.iter()) {
        *r = *limb;
    }

    E::Fr::from_repr(alpha_inv_repr).ok()
}

pub(crate) fn init_old_poseidon_params() -> poseidon_hash::bn256::Bn256PoseidonParams {
    poseidon_hash::bn256::Bn256PoseidonParams::new_checked_2_into_1()
}
