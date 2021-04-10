use franklin_crypto::{
    bellman::{
        plonk::better_better_cs::cs::{ConstraintSystem, PlonkConstraintSystemParams},
        Engine,
    },
    bellman::{Field, PrimeField, SynthesisError},
    plonk::circuit::allocated_num::AllocatedNum,
    plonk::circuit::{
        allocated_num::Num, custom_rescue_gate::apply_5th_power,
        linear_combination::LinearCombination,
    },
};

// Substitution box is non-linear part of permutation function.
// It basically computes 5th power of each element in the state.
// Poseidon uses partial sbox which basically computes power of
// single element of state. If constraint system has support of
// custom gate then computation costs only single gate.
pub(crate) fn sbox_quintic<E: Engine, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    prev_state: &mut [LinearCombination<E>],
) -> Result<(), SynthesisError> {
    let state_as_nums: Vec<Result<Num<E>, SynthesisError>> = prev_state
        .iter_mut()
        .map(|s| s.clone().into_num(cs))
        .collect();
    match CS::Params::HAS_CUSTOM_GATES == true && CS::Params::STATE_WIDTH >= 4 {
        true => {
            for (s, s_num) in prev_state.iter_mut().zip(state_as_nums) {
                *s = match s_num? {
                    Num::Variable(var) => LinearCombination::from(apply_5th_power(cs, &var, None)?),
                    Num::Constant(c) => {
                        let tmp = c.pow([5u64]);
                        let mut lc = LinearCombination::zero();
                        lc.add_assign_constant(tmp);
                        lc
                    }
                };
            }
            Ok(())
        }
        false => {
            for (s, s_num) in prev_state.iter_mut().zip(state_as_nums) {
                let s_num = s_num?.get_variable();
                let squa = s_num.square(cs)?;
                let quad = squa.square(cs)?;
                let quin = quad.mul(cs, &s_num)?;
                *s = LinearCombination::from(quin);
            }
            Ok(())
        }
    }
}
// This function computes power of inverse of alpha to each element of state.
// By custom gate support, it costs only single gate. Under the hood, it proves
// that 5th power of each element of state is equal to itself.(x^(1/5)^5==x)
pub(crate) fn sbox_quintic_inv<E: Engine, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    alpha_inv: E::Fr,
    prev_state: &mut [LinearCombination<E>],
) -> Result<(), SynthesisError> {
    let state_as_nums: Vec<Result<Num<E>, SynthesisError>> = prev_state
        .iter_mut()
        .map(|s| s.clone().into_num(cs))
        .collect();
    match CS::Params::HAS_CUSTOM_GATES == true && CS::Params::STATE_WIDTH >= 4 {
        false => unimplemented!(),
        true => {
            for (s, s_num) in prev_state.iter_mut().zip(state_as_nums) {
                // x^(1/alpha)
                let alpha_inv_alloc = AllocatedNum::alloc(cs, || {
                    let original = s_num?.get_value().expect("");
                    let result = original.pow(alpha_inv.into_repr());
                    Ok(result)
                })?;
                let _ = apply_5th_power(cs, &alpha_inv_alloc, None)?;
                *s = LinearCombination::from(alpha_inv_alloc);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use franklin_crypto::{
        bellman::bn256::{Bn256, Fr},
        bellman::PrimeField,
        plonk::circuit::{allocated_num::AllocatedNum, linear_combination::LinearCombination},
    };
    use rand::Rand;

    use crate::tests::{init_cs, init_rng};

    use super::*;

    #[test]
    fn test_sbox_quintic_with_custom_gate() {
        let cs = &mut init_cs::<Bn256>();
        let rng = &mut init_rng();

        let a = Fr::rand(rng);
        let a_num = AllocatedNum::alloc(cs, || Ok(a)).expect("valid el");
        let a_lc = LinearCombination::from(a_num);

        let n = 3;
        for _ in 0..n {
            let _ = sbox_quintic(cs, &mut [a_lc.clone()]).expect("5th apply successfu");
        }
        // cs.finalize();
        // assert!(cs.is_satisfied());

        // println!(
        //     "quintic sbox takes {} gates with custom gate for {} iteration ",
        //     cs.n(),
        //     n
        // );
    }
    #[test]
    fn test_sbox_alpha_without_custom_gate() {
        // TODO: use CS which has no custom gate support
        let cs = &mut init_cs::<Bn256>();
        let rng = &mut init_rng();

        let a = Fr::rand(rng);
        let a_num = AllocatedNum::alloc(cs, || Ok(a)).unwrap();
        let b = Fr::rand(rng);
        let b_num = AllocatedNum::alloc(cs, || Ok(b)).unwrap();
        let b_lc = LinearCombination::from(b_num);
        let mut a_lc = LinearCombination::from(a_num);
        a_lc.add_assign(&b_lc);
        let _alpha = Fr::from_str("5").unwrap();

        let n = 2;
        for _ in 0..n {
            let _ = sbox_quintic(cs, &mut [a_lc.clone()]).unwrap();
        }
        // cs.finalize();
        // assert!(cs.is_satisfied());

        // println!(
        //     "quintic sbox takes {} gates without custom gate for {} iteration ",
        //     cs.n(),
        //     n
        // );
    }
    #[test]
    fn test_sbox_alpha_inv_without_custom_gate() {
        let cs = &mut init_cs::<Bn256>();
        let rng = &mut init_rng();

        let a = Fr::rand(rng);
        let a_num = AllocatedNum::alloc(cs, || Ok(a)).unwrap();
        let a_lc = LinearCombination::from(a_num);
        let _alpha = Fr::from_str("5").unwrap();
        let alpha_inv = crate::common::utils::compute_gcd::<Bn256>(5u64);
        let _ = sbox_quintic_inv(cs, alpha_inv.expect("inverse of alpha"), &mut [a_lc]).unwrap();

        // cs.finalize();
        // assert!(cs.is_satisfied());

        // println!("power sbox takes {} gates without custom gate", cs.n());
    }
}
