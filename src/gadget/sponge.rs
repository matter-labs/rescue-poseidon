use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem,
    plonk::circuit::linear_combination::LinearCombination,
};
use franklin_crypto::{
    bellman::{Engine, Field, SynthesisError},
    plonk::circuit::{allocated_num::Num, boolean::Boolean},
};

// This traits will be moved into franklin, thats why we use different trait.
pub trait GadgetSpongeParams {
    fn rate(&self) -> usize;
}
pub trait GadgetSpongeState<E: Engine> {
    fn state_as_ref(&self) -> &[LinearCombination<E>];
    fn storage_as_ref(&self) -> &[Num<E>];
    fn state_as_mut(&mut self) -> &mut [LinearCombination<E>];
    fn storage_as_mut(&mut self) -> &mut Vec<Num<E>>;
}

pub trait GadgetSpongePermutation<E: Engine> {
    fn permutation<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        should_permute: &Boolean,
    ) -> Result<(), SynthesisError>;
}

#[derive(Clone, Debug)]
pub enum SpongeModes<E: Engine> {
    // Standard mode is stateless
    Standard,
    // Duplex is statefull and maximum number of element "l" one can request
    // is equal to rate parameter.
    Duplex(Vec<Num<E>>),
}

pub trait GadgetSpongeMode<E: Engine> {
    fn get_mode(&self) -> SpongeModes<E>;
    fn update_mode(&mut self, mode: SpongeModes<E>);
}

pub trait StatefulSpongeGadget<E: Engine>:
    GadgetSpongeState<E>
    + GadgetSpongePermutation<E>
    + GadgetSpongeParams
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
        input: Num<E>,
    ) -> Result<(), SynthesisError> {
        let rate = self.rate();
        if self.storage_as_ref().len() < rate {
            self.storage_as_mut().push(input);
        }

        match self.get_mode() {
            SpongeModes::Standard =>  {
                if self.storage_as_ref().len() == rate {
                    permute::<_, CS, _>(cs, self)
                }else{
                    Ok(())
                }
            },
            SpongeModes::Duplex(_) => {
                // If state already squeezed then discard buffer. We don't need to
                // accumulate any value here because we alread stored in top of function
                self.update_mode(SpongeModes::Duplex(Vec::with_capacity(self.rate())));
                Ok(())
            }
        }
    }

    fn absorb_multi<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        input: &[Num<E>],
    ) -> Result<(), SynthesisError> {
        let rate = self.rate();
        assert!(!input.is_empty());
        assert_eq!(input.len() % rate, 0);
        for i in 0..input.len() / rate {
            for value in input.iter().skip(i * rate).take(rate) {
                self.absorb(cs, *value)?;
            }
        }

        Ok(())
    }

    fn squeeze<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
    ) -> Result<Vec<Num<E>>, SynthesisError> {
        assert!(
            self.storage_as_ref().is_empty(),
            "storage should be empty which also means permutation happened"
        );
        assert!(
            self.state_as_ref()
                .iter()
                .all(|el| !el.get_value().expect("").is_zero()),
            "state elements should not equal to zero"
        );

        let output = self.state_as_ref()[..self.rate()]
            .iter()
            .map(|s| s.clone().into_num(cs).expect("should get a num"))
            .collect();
        Ok(output)
    }

    fn squeeze_single<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
    ) -> Result<Num<E>, SynthesisError> {
        match self.get_mode() {
            SpongeModes::Standard => {
                unimplemented!("Only duplex sponge can squeeze single element")
            }
            SpongeModes::Duplex(mut buffer) => {
                // use already squeezed values
                if buffer.len() > 0 {
                    let out = buffer.remove(0);
                    self.update_mode(SpongeModes::Duplex(buffer));
                    return Ok(out);
                }
                // at least one element should have been absorbed
                assert!(self.storage_as_ref().len() >= 1);
                // pad elements in order to run a permutation
                // TODO:  PaddingStrategy enum can be used here
                while self.storage_as_ref().len() % self.rate() != 0 {
                    self.storage_as_mut().push(Num::Constant(E::Fr::one()));
                }
                permute(cs, self)?;
                // squeeze from state
                let mut output: Vec<Num<E>> = self.state_as_ref()[..self.rate()]
                    .iter()
                    .map(|s| s.clone().into_num(cs).expect("should get a num"))
                    .collect();
                // output single element
                let out = output.remove(0);
                self.update_mode(SpongeModes::Duplex(output));
                Ok(out)
            }
        }
    }

    fn squeeze_multi<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        len: usize,
    ) -> Result<Vec<Num<E>>, SynthesisError> {
        assert!(self.storage_as_ref().is_empty(), "storage should be empty");
        assert!(self
            .state_as_ref()
            .iter()
            .all(|el| !el.get_value().expect("").is_zero()),);
        assert!(
            len > 2,
            "length of requested output should be greater than rate param"
        );
        let mut output = vec![];
        let should_permute = Boolean::Constant(true);
        while output.len() < len {
            for i in 0..self.rate() {
                output.push(self.state_as_ref()[i].clone().into_num(cs)?);
            }

            if output.len() < len {
                self.permutation(cs, &should_permute)?;
                self.storage_as_mut().truncate(0);
            }
        }
        assert!(output.len() == len);

        Ok(output)
    }

    fn reset(&mut self) {
        self.storage_as_mut().truncate(0);
        self.state_as_mut()
            .iter_mut()
            .for_each(|s| *s = LinearCombination::zero());
    }
}

fn permute<E: Engine, CS: ConstraintSystem<E>, S: StatefulSpongeGadget<E>>(
    cs: &mut CS,
    sponge: &mut S,
) -> Result<(), SynthesisError> {
    let storage_values = sponge.storage_as_ref().to_vec();
    for (value, state) in storage_values.iter().zip(sponge.state_as_mut().iter_mut()) {
        state.add_assign_number_with_coeff(value, E::Fr::one());
    }
    sponge.permutation(cs, &Boolean::Constant(true))?;
    sponge.storage_as_mut().truncate(0);
    Ok(())
}

#[macro_export]
macro_rules! sponge_gadget_impl {
    ($hasher_name:ty) => {
        impl<E: Engine> StatefulSpongeGadget<E> for $hasher_name {}

        impl<E: Engine> GadgetSpongeParams for $hasher_name {
            fn rate(&self) -> usize {
                self.params.rate
            }
        }

        impl<E: Engine> GadgetSpongeState<E> for $hasher_name {
            fn state_as_ref(&self) -> &[LinearCombination<E>] {
                self.state.as_ref()
            }
            fn storage_as_ref(&self) -> &[Num<E>] {
                self.tmp_storage.as_ref()
            }
            fn state_as_mut(&mut self) -> &mut [LinearCombination<E>] {
                self.state.as_mut()
            }
            fn storage_as_mut(&mut self) -> &mut Vec<Num<E>> {
                self.tmp_storage.as_mut()
            }
        }

        impl<E: Engine> GadgetSpongeMode<E> for $hasher_name {
            fn get_mode(&self) -> SpongeModes<E> {
                self.sponge_mode.to_owned()
            }
            fn update_mode(&mut self, mode: SpongeModes<E>) {
                self.sponge_mode = mode;
            }
        }
    };
}
