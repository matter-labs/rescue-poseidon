use franklin_crypto::{
    bellman::{
        plonk::better_better_cs::cs::{ArithmeticTerm, ConstraintSystem, MainGateTerm, PlonkConstraintSystemParams},
        Engine,
    },
    bellman::{Field, SynthesisError},
    plonk::circuit::allocated_num::AllocatedNum,
    plonk::circuit::{
        allocated_num::Num, custom_rescue_gate::apply_5th_power,
        linear_combination::LinearCombination,
    },
};

use crate::traits::Sbox;

// Substitution box is non-linear part of permutation function.
// It basically computes 5th power of each element in the state.
// Poseidon uses partial sbox which basically computes power of
// single element of state. If constraint system has support of
// custom gate then computation costs only single gate.
// TODO use const generics here
pub(crate) fn sbox<E: Engine, CS: ConstraintSystem<E>, const WIDTH: usize>(
    cs: &mut CS,
    power: &Sbox,
    prev_state: &mut [LinearCombination<E>; WIDTH],
    state_range: Option<std::ops::Range<usize>>,
    use_custom_gate: bool,
) -> Result<(), SynthesisError> {
    let use_custom_gate =
        use_custom_gate && CS::Params::HAS_CUSTOM_GATES == true && CS::Params::STATE_WIDTH >= 4;
    match power {
        Sbox::Alpha(alpha) => sbox_alpha(
            cs,
            alpha,
            prev_state,
            state_range.expect("full state not partial"),
            use_custom_gate,
        ),
        Sbox::AlphaInverse(alpha_inv) => {
            // TODO
            // assert!(
            //     state_range.is_none(),
            //     "partial sbox doesn't supported in inverse direction"
            // );
            sbox_alpha_inv(cs, alpha_inv, prev_state, use_custom_gate)
        }
    }
}

fn sbox_alpha<E: Engine, CS: ConstraintSystem<E>, const WIDTH: usize>(
    cs: &mut CS,
    alpha: &u64,
    prev_state: &mut [LinearCombination<E>; WIDTH],
    state_range: std::ops::Range<usize>,
    use_custom_gate: bool,
) -> Result<(), SynthesisError> {
    if *alpha != 5u64 {
        unimplemented!("only 5th power is supported!")
    }
    for lc in prev_state[state_range].iter_mut() {
        match lc.clone().into_num(cs)? {
            Num::Constant(value) => {
                let result = value.pow(&[*alpha]);
                *lc = LinearCombination::zero();
                lc.add_assign_constant(result);
            }
            Num::Variable(ref value) => {
                let result = if use_custom_gate {
                    apply_5th_power(cs, value, None)?
                } else {
                    let square = value.square(cs)?;
                    let quad = square.square(cs)?;
                    quad.mul(cs, value)?
                };
                *lc = LinearCombination::from(result);
            }
        }
    }

    return Ok(());
}
// This function computes power of inverse of alpha to each element of state.
// By custom gate support, it costs only single gate. Under the hood, it proves
// that 5th power of each element of state is equal to itself.(x^(1/5)^5==x)
fn sbox_alpha_inv<E: Engine, CS: ConstraintSystem<E>, const WIDTH: usize, const N: usize>(
    cs: &mut CS,
    alpha_inv: &[u64; N],
    prev_state: &mut [LinearCombination<E>; WIDTH],
    use_custom_gate: bool,
) -> Result<(), SynthesisError> {
    for lc in prev_state.iter_mut() {
        match lc.clone().into_num(cs)? {
            Num::Constant(value) => {
                let result = value.pow(alpha_inv);
                *lc = LinearCombination::zero();
                lc.add_assign_constant(result);
            }
            Num::Variable(ref value) => {
                let powered = AllocatedNum::alloc(cs, || {
                    let base = value.get_value().expect("value");
                    
                    let result = base.pow(alpha_inv);
                    Ok(result)
                })?;

                if use_custom_gate {
                    let _  = apply_5th_power(cs, &powered, Some(*value))?; 
                } else {
                    let squared = powered.square(cs)?;
                    let quad = squared.square(cs)?;

                    let mut term = MainGateTerm::<E>::new();
                    let fifth_term = ArithmeticTerm::from_variable(quad.get_variable())
                        .mul_by_variable(powered.get_variable());
                    let el_term = ArithmeticTerm::from_variable(value.get_variable());
                    term.add_assign(fifth_term);
                    term.sub_assign(el_term);
                    cs.allocate_main_gate(term)?;
                };
                *lc = LinearCombination::from(powered);
            }
        }
    }

    return Ok(());
}

#[cfg(test)]
mod test {
    use franklin_crypto::{
        bellman::{bn256::Bn256, PrimeField},
        plonk::circuit::linear_combination::LinearCombination,
    };
    use std::convert::TryInto;

    use super::super::tests::test_inputs;
    use crate::tests::{init_cs, init_cs_no_custom_gate};

    use super::*;

    fn test_sbox<E: Engine, CS: ConstraintSystem<E>, const N: usize>(
        cs: &mut CS,
        power: Sbox,
        number_of_rounds: usize,
        use_custom_gate: bool,
        use_partial_state: bool,
        use_allocated: bool,
    ) {
        let (mut state, state_as_nums) = test_inputs::<E, _, N>(cs, use_allocated);

        let mut state_as_lc = vec![];
        for num in std::array::IntoIter::new(state_as_nums) {
            let lc = LinearCombination::from(num);
            state_as_lc.push(lc);
        }
        let mut state_as_lc: [LinearCombination<E>; N] = state_as_lc.try_into().expect("array");

        let state_range = if use_partial_state {
            Some(0..1)
        } else {
            Some(0..N)
        };

        assert_eq!(state_range.as_ref().unwrap().len(), N);

        for _ in 0..number_of_rounds {
            crate::common::sbox::sbox::<E>(&power, &mut state);
            sbox(
                cs,
                &power,
                &mut state_as_lc,
                state_range.clone(),
                use_custom_gate,
            )
            .expect("5th apply successfu");
        }
    }
    #[test]
    fn test_sbox_quintic() {
        let cs = &mut init_cs::<Bn256>();

        const INPUT_LENGTH: usize = 3;
        const NUM_ROUNDS: usize = 1;

        let alpha = Sbox::Alpha(5);

        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha.clone(), NUM_ROUNDS, true, false, true); // variable inputs + custom gate
        println!(
            "quintic sbox takes {} gates with custom gate and variable inputs for {} iteration ",
            cs.n(),
            NUM_ROUNDS,
        );
        let mut end = cs.n();
        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha.clone(), NUM_ROUNDS, true, false, false); // constant inputs + custom gate
        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha.clone(), NUM_ROUNDS, false, false, true); // variable inputs + no custom gate
        println!(
            "quintic sbox takes {} gates without custom gate and variable inputs for {} iteration ",
            cs.n() - end,
            NUM_ROUNDS,
        );
        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha, NUM_ROUNDS, false, false, false); // constant inputs + no custom gate        

        cs.finalize();
        assert!(cs.is_satisfied());
    }

    #[test]
    fn test_sbox_quintic_inv() {
        let cs = &mut init_cs::<Bn256>();

        const INPUT_LENGTH: usize = 3;
        const NUM_ROUNDS: usize = 1;
        let alpha = 5;
        let alpha_inv = Sbox::AlphaInverse(compute_inverse_alpha::<Bn256, 4>(alpha));

        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha_inv.clone(), NUM_ROUNDS, true, false, true); // variable inputs + custom gate
        println!(
            "quintic inverse sbox takes {} gates with custom gate and variable inputs for {} iteration ",
            cs.n(),
            NUM_ROUNDS,
        );
        let mut end = cs.n();
        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha_inv.clone(), NUM_ROUNDS, true, false, false); // constant inputs + custom gate        
        end = cs.n();
        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha_inv.clone(), NUM_ROUNDS, false, false, true); // variable inputs + no custom gate
        println!(
            "quintic inverse sbox takes {} gates without custom gate and variable inputs for {} iteration ",
            cs.n() - end,
            NUM_ROUNDS,
        );
        end = cs.n();
        test_sbox::<_, _, INPUT_LENGTH>(cs, alpha_inv, NUM_ROUNDS, false, false, false); // constant inputs + no custom gate
        cs.finalize();
        assert!(cs.is_satisfied());
    }

    fn compute_inverse_alpha<E: Engine, const N: usize>(alpha: u64) -> [u64; N] {
        crate::common::utils::compute_gcd::<E, N>(alpha).expect("inverse of alpha")        
    }

}
