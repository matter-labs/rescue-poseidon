use crate::{
    common::domain_strategy::DomainStrategy,
    traits::{HashFamily, HashParams},
};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::allocated_num::Num,
};
use franklin_crypto::{bellman::Field, plonk::circuit::boolean::Boolean};
use franklin_crypto::{
    bellman::{Engine, SynthesisError},
    plonk::circuit::linear_combination::LinearCombination,
};
use std::convert::TryInto;

pub fn circuit_generic_hash<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
    const LENGTH: usize,
>(
    cs: &mut CS,
    input: &[Num<E>; LENGTH],
    params: &P,
) -> Result<[Num<E>; RATE], SynthesisError> {
    CircuitGenericSponge::hash(cs, input, params)
}

#[derive(Clone)]
enum SpongeMode<E: Engine, const RATE: usize> {
    Absorb([Option<Num<E>>; RATE]),
    Squeeze([Option<Num<E>>; RATE]),
}

#[derive(Clone)]
pub struct CircuitGenericSponge<E: Engine, const RATE: usize, const WIDTH: usize> {
    state: [LinearCombination<E>; WIDTH],
    mode: SpongeMode<E, RATE>,
}

impl<'a, E: Engine, const RATE: usize, const WIDTH: usize> CircuitGenericSponge<E, RATE, WIDTH> {
    pub fn new() -> Self {
        let state = (0..WIDTH)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<_>>()
            .try_into()
            .expect("constant array");

        Self {
            state,
            mode: SpongeMode::Absorb([None; RATE]),
        }
    }

    pub fn hash<CS: ConstraintSystem<E>, P: HashParams<E, RATE, WIDTH>>(
        cs: &mut CS,
        input: &[Num<E>],
        params: &P,
    ) -> Result<[Num<E>; RATE], SynthesisError> {
        // init state
        let mut state: [LinearCombination<E>; WIDTH] = (0..WIDTH)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("constant array of LCs");

        let domain_strategy = DomainStrategy::CustomFixedLength;
        // specialize capacity
        let capacity_value = domain_strategy
            .compute_capacity::<E>(input.len(), RATE)
            .unwrap_or(E::Fr::zero());
        state
            .last_mut()
            .expect("last element")
            .add_assign_constant(capacity_value);

        // compute padding values
        let padding_values = domain_strategy
            .generate_padding_values::<E>(input.len(), RATE)
            .iter()
            .map(|el| Num::Constant(*el))
            .collect::<Vec<Num<E>>>();

        // chain all values
        let mut padded_input = vec![];
        padded_input.extend_from_slice(input);
        padded_input.extend_from_slice(&padding_values);

        assert!(padded_input.len() % RATE == 0);

        // process each chunk of input
        for values in padded_input.chunks_exact(RATE) {
            absorb(
                cs,
                &mut state,
                values.try_into().expect("constant array"),
                params,
            )?;
        }
        // prepare output
        let mut output = [Num::Constant(E::Fr::zero()); RATE];
        for (o, s) in output.iter_mut().zip(state.iter()) {
            *o = s.clone().into_num(cs)?;
        }

        Ok(output)
    }

    pub fn absorb_multiple<CS: ConstraintSystem<E>, P: HashParams<E, RATE, WIDTH>>(
        &mut self,
        cs: &mut CS,
        input: &[Num<E>],
        params: &P,
    ) -> Result<(), SynthesisError> {
        for inp in input.into_iter() {
            self.absorb(cs, *inp, params)?
        }

        Ok(())
    }

    pub fn absorb<CS: ConstraintSystem<E>, P: HashParams<E, RATE, WIDTH>>(
        &mut self,
        cs: &mut CS,
        input: Num<E>,
        params: &P,
    ) -> Result<(), SynthesisError> {
        match self.mode {
            SpongeMode::Absorb(ref mut buf) => {
                // push value into buffer
                for el in buf.iter_mut() {
                    if el.is_none() {
                        // we still have empty room for values
                        *el = Some(input);
                        return Ok(());
                    }
                }

                // buffer is filled so unwrap them
                let mut unwrapped_buffer = [Num::Constant(E::Fr::zero()); RATE];
                for (a, b) in unwrapped_buffer.iter_mut().zip(buf.iter_mut()) {
                    if let Some(val) = b {
                        *a = *val;
                        *b = None; // kind of resetting buffer
                    }
                }

                // here we can absorb values. run round function implicitly there
                absorb::<_, _, P, RATE, WIDTH>(cs, &mut self.state, &mut unwrapped_buffer, params)?;

                // absorb value
                buf[0] = Some(input);
            }
            SpongeMode::Squeeze(_) => {
                // we don't need squeezed values so switching to absorbing mode is fine
                let mut buf = [None; RATE];
                buf[0] = Some(input);
                self.mode = SpongeMode::Absorb(buf)
            }
        }

        Ok(())
    }

    /// Apply padding manually especially when single absorb called single/many times
    pub fn pad_if_necessary(&mut self) {
        match self.mode {
            SpongeMode::Absorb(ref mut buf) => {
                let unwrapped_buffer_len = buf.iter().filter(|el| el.is_some()).count();
                // compute padding values
                let padding_strategy = DomainStrategy::CustomVariableLength;
                let padding_values =
                    padding_strategy.generate_padding_values::<E>(unwrapped_buffer_len, RATE);
                let mut padding_values_it = padding_values.iter().cloned();

                for b in buf {
                    if b.is_none() {
                        *b = Some(Num::Constant(padding_values_it.next().expect("next elm")))
                    }
                }
                assert!(padding_values_it.next().is_none());
            }
            SpongeMode::Squeeze(_) => (),
        }
    }
    
    pub fn squeeze<CS: ConstraintSystem<E>, P: HashParams<E, RATE, WIDTH>>(
        &mut self,
        cs: &mut CS,
        params: &P,
    ) -> Result<Option<Num<E>>, SynthesisError> {
        loop {
            match self.mode {
                SpongeMode::Absorb(ref mut buf) => {
                    // buffer may not be filled fully so we may need padding.
                    let mut unwrapped_buffer = vec![];
                    for el in buf {
                        if let Some(value) = el {
                            unwrapped_buffer.push(*value);
                        }
                    }

                    if unwrapped_buffer.len() != RATE {
                        // processing buffer was done and we need padding                        
                        return Ok(None);
                    }

                    // make input array
                    let mut all_inputs = [Num::Constant(E::Fr::zero()); RATE];
                    for (a, b) in all_inputs.iter_mut().zip(unwrapped_buffer) {
                        *a = b;
                    }

                    // permute state
                    absorb(cs, &mut self.state, &all_inputs, params)?;

                    // push values into squeezing buffer for later squeezing
                    let mut squeeze_buffer = [None; RATE];
                    for (s, b) in self.state[..RATE].iter().zip(squeeze_buffer.iter_mut()) {
                        *b = Some(s.clone().into_num(cs)?)
                    }

                    // we are switching squeezing mode so we can ignore to reset absorbing buffer
                    self.mode = SpongeMode::Squeeze(squeeze_buffer);
                }
                SpongeMode::Squeeze(ref mut buf) => {
                    for el in buf {
                        if let Some(value) = el.take() {
                            return Ok(Some(value));
                        }
                    }
                    return Ok(None);
                }
            };
        }
    }
}

fn absorb<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    state: &mut [LinearCombination<E>; WIDTH],
    input: &[Num<E>; RATE],
    params: &P,
) -> Result<(), SynthesisError> {
    for (v, s) in input.iter().zip(state.iter_mut()) {
        s.add_assign_number_with_coeff(v, E::Fr::one());
    }
    circuit_generic_round_function(cs, state, params)
}

pub fn circuit_generic_round_function<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    state: &mut [LinearCombination<E>; WIDTH],
    params: &P,
) -> Result<(), SynthesisError> {
    match params.hash_family() {
        HashFamily::Rescue => super::rescue::circuit_rescue_round_function(cs, params, state),
        HashFamily::Poseidon => super::poseidon::circuit_poseidon_round_function(cs, params, state),
        HashFamily::RescuePrime => {
            super::rescue_prime::gadget_rescue_prime_round_function(cs, params, state)
        }
    }
}

pub fn circuit_generic_round_function_conditional<
    E: Engine,
    CS: ConstraintSystem<E>,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    state: &mut [LinearCombination<E>; WIDTH],
    execute: &Boolean,
    params: &P,
) -> Result<(), SynthesisError> {
    match params.hash_family() {
        HashFamily::Rescue => super::rescue::circuit_rescue_round_function(cs, params, state),
        HashFamily::Poseidon => super::poseidon::circuit_poseidon_round_function(cs, params, state),
        HashFamily::RescuePrime => {
            super::rescue_prime::gadget_rescue_prime_round_function(cs, params, state)
        }
    }
}
