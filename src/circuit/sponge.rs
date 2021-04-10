use std::convert::TryInto;

use super::traits::SpongeGadget;
use crate::traits::HashFamily;
use crate::{sponge::SpongeModes, traits::HashParams};
use franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use franklin_crypto::{
    bellman::{Engine, Field, SynthesisError},
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};

pub struct GenericSpongeGadget<
    'a,
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    const STATE_WIDTH: usize,
    const RATE: usize,
> {
    params: &'a P,
    state: [LinearCombination<E>; STATE_WIDTH],
    mode: SpongeModes,
}

impl<
        'a,
        E: Engine,
        P: HashParams<E, STATE_WIDTH, RATE>,
        const STATE_WIDTH: usize,
        const RATE: usize,
    > From<&'a P> for GenericSpongeGadget<'a, E, P, STATE_WIDTH, RATE>
{
    fn from(params: &'a P) -> Self {
        let state = (0..STATE_WIDTH)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("constant array of LCs");
        Self {
            params,
            state,
            mode: SpongeModes::Standard(false),
        }
    }
}

impl<
        'a,
        E: Engine,
        P: HashParams<E, STATE_WIDTH, RATE>,
        const STATE_WIDTH: usize,
        const RATE: usize,
    > SpongeGadget<E, STATE_WIDTH, RATE> for GenericSpongeGadget<'a, E, P, STATE_WIDTH, RATE>
{
    fn specialize(
        &mut self,
        capacity_value: Option<LinearCombination<E>>, // TODO Fr?
    ) -> Result<(), SynthesisError> {
        let value = capacity_value.unwrap_or(LinearCombination::zero());
        if let Some(last_el) = self.state.last_mut() {
            *last_el = value
        }

        Ok(())
    }

    fn absorb<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        input: &[Num<E>],
    ) -> Result<(), SynthesisError> {
        assert!(!input.is_empty());

        match self.mode {
            SpongeModes::Standard(ref mut is_absorbed) => {
                assert_eq!(
                    input.len() % RATE,
                    0,
                    "input length is not multiple of rate"
                );
                assert!(!*is_absorbed, "Sponge should be in in absorbtion phase");
                for elems in input.chunks_exact(RATE) {
                    for (value, state) in elems.iter().zip(self.state.iter_mut()) {
                        state.add_assign_number_with_coeff(value, E::Fr::one());
                    }
                    generic_round_function_gadget(cs, self.params, &mut self.state)?;
                    *is_absorbed = true;
                }
                
            }
            SpongeModes::Duplex(ref mut is_absorbed) => {
                assert!(!*is_absorbed, "Sponge should be in in absorbtion phase");
                assert!(
                    input.len() <= RATE,
                    "duplex sponge can absorb max rate elems"
                );
                for (value, state) in input.iter().zip(self.state.iter_mut()) {
                    state.add_assign_number_with_coeff(value, E::Fr::one());
                }
                generic_round_function_gadget(cs, self.params, &mut self.state)?;
                *is_absorbed = true;
            }
        }
        Ok(())
    }

    fn squeeze<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        number_of_elems: Option<usize>,
    ) -> Result<Vec<Num<E>>, SynthesisError> {
        let mut out = vec![];

        match self.mode {
            SpongeModes::Standard(ref mut is_absorbed) => {
                assert!(*is_absorbed, "Sponge should be in in squeezing phase");
                if let Some(number_of_elems) = number_of_elems {
                    if number_of_elems <= RATE {
                        out.extend_from_slice(&self.state[..RATE]);
                    } else {
                        let original_number_of_elems = number_of_elems;

                        let number_of_iters = if number_of_elems % RATE != 0 {
                            (number_of_elems + (RATE - (number_of_elems % RATE))) / RATE
                        } else {
                            number_of_elems / RATE
                        };

                        for _ in 0..number_of_iters {
                            out.extend_from_slice(&self.state[..RATE]);
                            generic_round_function_gadget(cs, self.params, &mut self.state)?;
                        }

                        out.truncate(original_number_of_elems);
                    }
                } else {
                    out.extend_from_slice(&self.state[..RATE]);
                }
                *is_absorbed = false;
                self.reset();
            }

            SpongeModes::Duplex(ref mut is_absorbed) => {
                assert!(*is_absorbed, "Sponge should be in in squeezing phase");
                let number_of_elems = if let Some(number_of_elems) = number_of_elems {
                    assert!(
                        number_of_elems <= RATE,
                        "duplex sponge squeeze only as much as rate parameter"
                    );
                    number_of_elems
                } else {
                    RATE
                };
                out.extend_from_slice(&self.state[..number_of_elems]);
                *is_absorbed = false;
            }
        }

        let out: Vec<Num<E>> = out
            .iter()
            .map(|s| s.clone().into_num(cs).expect("a num"))
            .collect();

        Ok(out)
    }

    fn reset(&mut self) {
        self.state
            .iter_mut()
            .for_each(|s| *s = LinearCombination::zero());
    }
}

pub fn generic_round_function_gadget<
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
    match params.hash_family() {
        HashFamily::Rescue => super::rescue::gadget_rescue_round_function(cs, params, state),
        HashFamily::Poseidon => super::poseidon::gadget_poseidon_round_function(cs, params, state),
        HashFamily::RescuePrime => {
            super::rescue_prime::gadget_rescue_prime_round_function(cs, params, state)
        }
    }
}
