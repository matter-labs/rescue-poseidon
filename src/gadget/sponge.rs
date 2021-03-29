use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem,
    plonk::circuit::linear_combination::LinearCombination,
};
use franklin_crypto::{
    bellman::{Engine, Field, SynthesisError},
    plonk::circuit::{allocated_num::Num, boolean::Boolean},
};

pub trait GadgetSpongeState<E: Engine, const S: usize> {
    fn state_as_ref(&self) -> &[LinearCombination<E>; S];
    fn state_as_mut(&mut self) -> &mut [LinearCombination<E>; S];
}

pub trait GadgetSpongePermutation<E: Engine> {
    fn permutation<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        should_permute: &Boolean,
    ) -> Result<(), SynthesisError>;
}

#[derive(Clone, Debug)]
pub enum SpongeModes{
    // Standard mode is stateless
    Standard(bool),
    // Duplex is statefull and maximum number of element "l" one can request
    // is equal to rate parameter.
    Duplex(bool),
}

pub trait GadgetSpongeMode<E: Engine> {
    fn get_mode(&self) -> SpongeModes;
    fn update_mode(&mut self, mode: SpongeModes);
}

pub trait StatefulSpongeGadget<E: Engine, const S: usize, const R: usize>:
    GadgetSpongeState<E, S>
    + GadgetSpongePermutation<E>    
    + GadgetSpongeMode<E>
    + Default
{
    fn specialize(&mut self, capacity_value: Option<LinearCombination<E>>) {
        let state = self.state_as_mut();
        let value = capacity_value.unwrap_or(LinearCombination::zero());
        if let Some(last_el) = state.last_mut() {
            *last_el = value
        }
    }

    fn absorb<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        input: &[Num<E>],
    ) -> Result<(), SynthesisError> {
        assert!(!input.is_empty());
        let rate = R;        
        

        match self.get_mode() {
            SpongeModes::Standard(is_absorbed) =>  {
                assert_eq!(
                    input.len() % rate,
                    0,
                    "input length is not multiple of rate"
                );
                assert!(!is_absorbed, "Sponge should be in in absorbtion phase");
                for elems in input.chunks_exact(rate) {
                    for (value, state) in elems.iter().zip(self.state_as_mut().iter_mut()) {
                        state.add_assign_number_with_coeff(value, E::Fr::one());
                    }
                    self.permutation(cs, &Boolean::constant(true));
                    self.update_mode(SpongeModes::Standard(true));
                }
            },
            SpongeModes::Duplex(is_absorbed) => {
                assert!(!is_absorbed, "Sponge should be in in absorbtion phase");
                assert!(
                    input.len() <= rate,
                    "duplex sponge can absorb max rate elems"
                );
                // If state already squeezed then discard buffer. We don't need to
                // accumulate any value here because we alread stored in top of function
                // TODO
                for (value, state) in input.iter().zip(self.state_as_mut().iter_mut()) {
                    state.add_assign_number_with_coeff(value, E::Fr::one());
                }
                self.permutation(cs, &Boolean::constant(true));
                self.update_mode(SpongeModes::Standard(true));

            }
        }
        Ok(())
    }


    fn squeeze<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        number_of_elems: Option<usize>
    ) -> Result<Vec<Num<E>>, SynthesisError> {
        let rate = R;

        let mut out = vec![];

        match self.get_mode() {
            SpongeModes::Standard(is_absorbed) => {
                assert!(is_absorbed, "Sponge should be in in squeezing phase");
                if let Some(number_of_elems) = number_of_elems {
                    if number_of_elems <= rate {
                        out.extend_from_slice(&self.state_as_ref()[..rate]);
                    } else {
                        let original_number_of_elems = number_of_elems;

                        let number_of_iters = if number_of_elems % rate != 0 {
                            (number_of_elems + (rate - (number_of_elems % rate))) / rate
                        } else {
                            number_of_elems / rate
                        };

                        for _ in 0..number_of_iters {
                            out.extend_from_slice(&self.state_as_ref()[..rate]);
                            self.permutation(cs, &Boolean::constant(true));
                        }

                        out.truncate(original_number_of_elems);
                    }
                } else {
                    out.extend_from_slice(&self.state_as_ref()[..rate]);
                }
                self.update_mode(SpongeModes::Standard(false));
                self.reset();
            }

            SpongeModes::Duplex(is_absorbed) => {
                assert!(is_absorbed, "Sponge should be in in squeezing phase");
                let number_of_elems = if let Some(number_of_elems) = number_of_elems {
                    assert!(
                        number_of_elems <= rate,
                        "duplex sponge squeeze only as much as rate parameter"
                    );
                    number_of_elems
                } else {
                    rate
                };

                out.extend_from_slice(&self.state_as_ref()[..number_of_elems]);
                self.update_mode(SpongeModes::Standard(false));
            }
        }

        let out: Vec<Num<E>> = out.iter().map(|s| s.clone().into_num(cs).expect("a num")).collect();

        Ok(out)
    }


    fn reset(&mut self) {
        self.state_as_mut()
            .iter_mut()
            .for_each(|s| *s = LinearCombination::zero());
    }
}


#[macro_export]
macro_rules! sponge_gadget_impl {
    ($hasher_name:ty) => {
        impl<E: Engine, const S: usize, const R: usize> StatefulSpongeGadget<E, S, R> for $hasher_name {}

        impl<E: Engine, const S: usize, const R: usize> GadgetSpongeState<E, S> for $hasher_name {
            fn state_as_ref(&self) -> &[LinearCombination<E>; S] {
                &self.state
            }
            fn state_as_mut(&mut self) -> &mut [LinearCombination<E>; S] {
                &mut self.state
            }
        }

        impl<E: Engine, const S: usize, const R: usize> GadgetSpongeMode<E> for $hasher_name {
            fn get_mode(&self) -> SpongeModes {
                self.sponge_mode.to_owned()
            }
            fn update_mode(&mut self, mode: SpongeModes) {
                self.sponge_mode = mode;
            }
        }
    };
}
