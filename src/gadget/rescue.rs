use super::sbox::*;
use super::{
    sponge::{
        GadgetSpongeMode, GadgetSpongePermutation, GadgetSpongeState,
        SpongeModes, StatefulSpongeGadget,
    },
    utils::matrix_vector_product,
};
use crate::{common::domain_strategy::DomainStrategy, sponge_gadget_impl, HasherParams};
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
    params: HasherParams<E, S, R>,
    _alpha: E::Fr,
    alpha_inv: E::Fr,
    sponge_mode: SpongeModes,
}

impl<E: Engine, const S: usize, const R: usize> Default for RescueGadget<E, S, R> {
    fn default() -> Self {
        let (params, _alpha, alpha_inv) = crate::rescue::params::rescue_params();
        let initial_state: [LinearCombination<E>; S] = (0..S)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("vector of lc");
        Self {
            state: initial_state,
            params,
            _alpha,
            alpha_inv: alpha_inv.expect("inverse of alpha"),
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}

sponge_gadget_impl!(RescueGadget<E, S, R>);

impl<E: Engine, const S: usize, const R: usize> GadgetSpongePermutation<E> for RescueGadget<E, S, R> {
    fn permutation<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        _should_permute: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.state
            .iter_mut()
            .zip(self.params.constants_of_round(0).iter())
            .for_each(|(s, c)| s.add_assign_constant(*c));

        for round in 0..2 * self.params.full_rounds {
            // apply sbox
            if round & 1 == 0 {
                sbox_quintic_inv::<E, _>(cs, self.alpha_inv, &mut self.state)?;
            } else {
                sbox_quintic(cs, &mut self.state)?;
            }
            // mds row
            self.state = matrix_vector_product(cs, &self.params.mds_matrix(), &self.state)?;

            // round constants
            for (s, c) in self
                .state
                .iter_mut()
                .zip(self.params.constants_of_round(round + 1).iter().cloned())
            {
                s.add_assign_constant(c);
            }
        }
        Ok(())
    }
}

impl<E: Engine, const S: usize, const R: usize> RescueGadget<E, S, R> {
    pub fn new() -> Self {
        let (params, _alpha, alpha_inv) = crate::rescue::params::rescue_params();
        let initial_state: [LinearCombination<E>; S] = (0..S)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("vector of lc");
        Self {
            state: initial_state,
            params,
            _alpha,
            alpha_inv: alpha_inv.expect("inverse of alpha"),
            sponge_mode: SpongeModes::Standard(false),
        }
    }
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

        let mut gadget = RescueGadget::<_, STATE_WIDTH, RATE>::new();
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
