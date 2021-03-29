use super::sponge::{
    GadgetSpongeMode, GadgetSpongePermutation, GadgetSpongeState, SpongeModes,
    StatefulSpongeGadget,
};
use super::{sbox::*, utils::matrix_vector_product};
use crate::sponge_gadget_impl;
use crate::{common::padding::PaddingStrategy, HasherParams};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::boolean::Boolean,
};

use franklin_crypto::bellman::SynthesisError;
use franklin_crypto::{
    bellman::Engine,
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};

use std::convert::TryInto;

pub fn rescue_prime_gadget_fixed_length<E, CS, const S: usize, const R: usize>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<Vec<Num<E>>, SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    super::hash::generic_hash::<E, _, RescuePrimeGadget<E, S, R>, S, R>(
        cs,
        input,
        PaddingStrategy::FixedLength,
    )
}

pub fn rescue_prime_gadget_var_length<E, CS, const S: usize, const R: usize>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<Vec<Num<E>>, SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    super::hash::generic_hash::<E, _, RescuePrimeGadget<E, S, R>, S, R>(
        cs,
        input,
        PaddingStrategy::VariableLength,
    )
}

pub fn rescue_prime_gadget<E, CS, const S: usize, const R: usize>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<Vec<Num<E>>, SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    super::hash::generic_hash::<E, _, RescuePrimeGadget<E, S, R>, S, R>(
        cs,
        input,
        PaddingStrategy::Custom,
    )
}

pub struct RescuePrimeGadget<E: Engine, const S: usize, const R: usize> {
    state: [LinearCombination<E>; S],
    params: HasherParams<E, S, R>,
    _alpha: E::Fr,
    alpha_inv: E::Fr,
    sponge_mode: SpongeModes,
}

impl<E: Engine, const S: usize, const R: usize> Default for RescuePrimeGadget<E, S, R> {
    fn default() -> Self {
        let (params, alpha, alpha_inv) = crate::rescue_prime::params::rescue_prime_params();
        let initial_state: [LinearCombination<E>; S] = (0..S)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("vector of lc");
        Self {
            state: initial_state,
            params,
            _alpha: alpha,
            alpha_inv,
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}

sponge_gadget_impl!(RescuePrimeGadget<E, S, R>);

impl<E: Engine, const S: usize, const R: usize> GadgetSpongePermutation<E> for RescuePrimeGadget<E, S, R> {
    // permutation happens in 9 rounds
    // first round is sparse and other 8 full rounds are full
    // total cost 2 + 3*2 + 8*3*(2+2) = 104
    fn permutation<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        _should_permute: &Boolean,
    ) -> Result<(), SynthesisError> {
        for round in 0..self.params.full_rounds {
            // apply sbox
            // each lc will have 3 terms but there will be 1 in first iteration
            // total cost 2 gate per state vars = 6
            sbox_quintic(cs, &mut self.state)?;

            // mul by mds
            self.state = matrix_vector_product(cs, &self.params.mds_matrix(), &self.state)?;

            // round constants
            let constants = self.params.constants_of_round(round);
            for (s, c) in self.state.iter_mut().zip(constants.iter().cloned()) {
                s.add_assign_constant(c);
            }
            // apply inverse sbox
            sbox_quintic_inv::<E, _>(cs, self.alpha_inv, &mut self.state)?;

            // mul by mds
            self.state = matrix_vector_product(cs, &self.params.mds_matrix(), &self.state)?;

            // round constants
            let constants = self.params.constants_of_round(round + 1);
            for (s, c) in self.state.iter_mut().zip(constants.iter().cloned()) {
                s.add_assign_constant(c);
            }
        }
        Ok(())
    }
}

impl<E: Engine, const S: usize, const R: usize> RescuePrimeGadget<E, S, R> {
    pub fn new() -> Self {
        let (params, alpha, alpha_inv) = crate::rescue_prime::params::rescue_prime_params();
        let initial_state: [LinearCombination<E>; S] = (0..S)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("vector of lc");
        Self {
            state: initial_state,            
            params,
            _alpha: alpha,
            alpha_inv,
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}

#[cfg(test)]
mod test {
    use franklin_crypto::{
        bellman::{plonk::better_better_cs::cs::ConstraintSystem, Field},
        plonk::circuit::allocated_num::AllocatedNum,
        plonk::circuit::allocated_num::Num,
    };

    use super::RescuePrimeGadget;
    use crate::sponge::StatefulSponge;
    use crate::tests::init_cs;
    use crate::{gadget::sponge::StatefulSpongeGadget, tests::init_rng};
    use franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
    use rand::Rand;

    #[test]
    fn test_rescue_prime_sponge_with_custom_gate() {
        const STATE_WIDTH: usize = 3;
        const RATE: usize = 2;
        let rng = &mut init_rng();
        let cs = &mut init_cs::<Bn256>();

        let mut inputs = vec![Fr::zero(); RATE];
        let mut inputs_as_num = vec![Num::Constant(Fr::zero()); 2];
        for (i1, i2) in inputs.iter_mut().zip(inputs_as_num.iter_mut()) {
            *i1 = Fr::rand(rng);
            *i2 = Num::Variable(AllocatedNum::alloc(cs, || Ok(*i1)).unwrap());
        }

        let mut rescue_prime_gadget = RescuePrimeGadget::<_, STATE_WIDTH, RATE>::new();
        rescue_prime_gadget
            .absorb(cs, &inputs_as_num)
            .unwrap();
        let gadget_output = rescue_prime_gadget.squeeze(cs, None).unwrap();
        // cs.finalize();
        // assert!(cs.is_satisfied());

        println!("number of gates {}", cs.n());
        println!("last step number {}", cs.get_current_step_number());

        // rescue prime original
        let mut rescue_prime = crate::rescue_prime::RescuePrimeHasher::<Bn256, STATE_WIDTH, RATE>::default();
        // TODO
        rescue_prime.absorb(&inputs);
        // let output = rescue_prime.squeeze();
        let output = rescue_prime.squeeze(None);


        for (sponge, gadget) in output.iter().zip(gadget_output.iter()) {
            assert_eq!(gadget.get_value().unwrap(), *sponge);
        }
    }
}
