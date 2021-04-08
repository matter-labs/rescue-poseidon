use super::sbox::*;
use super::{
    sponge::{
        GadgetSpongeMode, GadgetSpongePermutation, GadgetSpongeState, SpongeModes,
        StatefulSpongeGadget,
    },
    utils::matrix_vector_product,
};
use crate::{
    common::domain_strategy::DomainStrategy,
    common::params::HasherParams,
    rescue::{HashParams, RescueParams},
    sponge_gadget_impl,
};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::boolean::Boolean,
};

use franklin_crypto::bellman::SynthesisError;
use franklin_crypto::{
    bellman::Engine,
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};

use std::convert::TryInto;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter. Uses state-width=3 and rate=2.
pub fn rescue_gadget<E, CS>(cs: &mut CS, input: &[Num<E>]) -> Result<[Num<E>; 2], SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    inner_rescue_gadget::<_, _>(cs, input, DomainStrategy::CustomFixedLength)
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule. Uses state-width=3 and rate=2.
pub fn rescue_gadget_var_length<E, CS>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<[Num<E>; 2], SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    inner_rescue_gadget::<_, _>(cs, input, DomainStrategy::CustomVariableLength)
}

fn inner_rescue_gadget<E, CS>(
    cs: &mut CS,
    input: &[Num<E>],
    domain: DomainStrategy<2>,
) -> Result<[Num<E>; 2], SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;

    let result =
        super::hash::generic_hash::<E, _, RescueGadget<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            cs, input, domain,
        )?;

    Ok(result.try_into().expect("fixed length array"))
}

pub struct RescueGadget<E: Engine, const S: usize, const R: usize> {
    state: [LinearCombination<E>; S],
    params: RescueParams<E, S, R>,
    sponge_mode: SpongeModes,
}

impl<E: Engine, const S: usize, const R: usize> Default for RescueGadget<E, S, R> {
    fn default() -> Self {
        let initial_state: [LinearCombination<E>; S] = (0..S)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("vector of lc");
        Self {
            state: initial_state,
            params: RescueParams::default(),
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}

sponge_gadget_impl!(RescueGadget<E, S, R>);

impl<E: Engine, const S: usize, const R: usize> GadgetSpongePermutation<E>
    for RescueGadget<E, S, R>
{
    fn permutation<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        _should_permute: &Boolean,
    ) -> Result<(), SynthesisError> {
        rescue_circuit_round_function(cs, &self.params, &mut self.state)
    }
}

pub fn rescue_circuit_round_function<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, STATE_WIDTH, RATE>,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    cs: &mut CS,
    params: &P,
    state: &mut [LinearCombination<E>; STATE_WIDTH],
) -> Result<(), SynthesisError> {
    state
        .iter_mut()
        .zip(params.constants_of_round(0).iter())
        .for_each(|(s, c)| s.add_assign_constant(*c));

    for round in 0..2 * params.number_of_full_rounds() {
        // apply sbox
        if round & 1 == 0 {
            sbox_quintic_inv::<E, _>(cs, params.alpha_inv(), state)?;
        } else {
            sbox_quintic(cs, state)?;
        }
        // mds row
        // TODO remove mut from mds
        *state = matrix_vector_product(cs, &params.mds_matrix(), state)?;

        // round constants
        for (s, c) in state
            .iter_mut()
            .zip(params.constants_of_round(round + 1).iter().cloned())
        {
            s.add_assign_constant(c);
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use franklin_crypto::{
        bellman::plonk::better_better_cs::cs::ConstraintSystem, bellman::Field,
        plonk::circuit::allocated_num::AllocatedNum, plonk::circuit::allocated_num::Num,
    };

    use crate::sponge::StatefulSponge;
    use crate::tests::init_cs;
    use crate::{gadget::sponge::StatefulSpongeGadget, tests::init_rng};
    use franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
    use rand::Rand;

    use super::RescueGadget;
    #[test]
    fn test_rescue_sponge_with_custom_gate() {
        const STATE_WIDTH: usize = 3;
        const RATE: usize = 2;
        let rng = &mut init_rng();
        let cs = &mut init_cs::<Bn256>();

        let mut inputs = vec![Fr::zero(); 2];
        let mut inputs_as_num = vec![Num::Constant(Fr::zero()); 2];
        for (i1, i2) in inputs.iter_mut().zip(inputs_as_num.iter_mut()) {
            *i1 = Fr::rand(rng);
            *i2 = Num::Variable(AllocatedNum::alloc(cs, || Ok(*i1)).unwrap());
        }

        let mut gadget = RescueGadget::<_, STATE_WIDTH, RATE>::default();
        gadget.absorb(cs, &inputs_as_num).unwrap();
        let gadget_output = gadget.squeeze(cs, None).unwrap();

        // cs.finalize();
        assert!(cs.is_satisfied());

        println!("number of gates {}", cs.n());
        println!("last step number {}", cs.get_current_step_number());

        // rescue original
        let mut rescue = crate::rescue::RescueHasher::<Bn256, STATE_WIDTH, RATE>::default();
        rescue.absorb(&inputs);
        let output = rescue.squeeze(None);

        for (sponge, gadget) in output.iter().zip(gadget_output.iter()) {
            assert_eq!(gadget.get_value().unwrap(), *sponge);
        }
    }
}
