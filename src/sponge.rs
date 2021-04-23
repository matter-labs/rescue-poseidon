use crate::{common::domain_strategy::DomainStrategy, traits::HashParams};
use franklin_crypto::bellman::Engine;
use franklin_crypto::bellman::Field;
use std::convert::TryInto;

pub fn generic_hash<
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
    const LENGTH: usize,
>(
    params: &P,
    input: &[E::Fr; LENGTH],
) -> [E::Fr; RATE] {
    GenericSponge::hash(input, params)
}

#[derive(Clone)]
enum SpongeMode<E: Engine, const RATE: usize> {
    Absorb([Option<E::Fr>; RATE]),
    Squeeze([Option<E::Fr>; RATE]),
}

#[derive(Clone)]
pub struct GenericSponge<
    'a,
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
> {
    state: [E::Fr; WIDTH],
    params: &'a P,
    mode: SpongeMode<E, RATE>,
}

impl<'a, E: Engine, P: HashParams<E, RATE, WIDTH>, const RATE: usize, const WIDTH: usize>
    GenericSponge<'a, E, P, RATE, WIDTH>
{
    pub fn new_from_params(params: &'a P) -> Self {
        Self::new_from_params_and_state([E::Fr::zero(); WIDTH], params)
    }

    pub fn new_from_params_and_state(state: [E::Fr; WIDTH], params: &'a P) -> Self {
        Self {
            state,
            params,
            mode: SpongeMode::Absorb([None; RATE]),
        }
    }

    pub fn hash(input: &[E::Fr], params: &P) -> [E::Fr; RATE] {
        // init state
        let mut state = [E::Fr::zero(); WIDTH];

        let domain_strategy = DomainStrategy::CustomFixedLength;
        // specialize capacity
        let capacity_value = domain_strategy
            .compute_capacity::<E>(input.len(), RATE)
            .unwrap_or(E::Fr::zero());
        *state.last_mut().expect("last element") = capacity_value;

        // compute padding values
        let padding_values = domain_strategy.generate_padding_values::<E>(input.len(), RATE);

        // chain all values
        let mut padded_input = vec![];
        padded_input.extend_from_slice(input);
        padded_input.extend_from_slice(&padding_values);

        assert!(padded_input.len() % RATE == 0);

        // process each chunk of input
        for values in padded_input.chunks_exact(RATE) {
            absorb::<E, _, RATE, WIDTH>(
                &mut state,
                &values.try_into().expect("constant array"),
                params,
            );
        }
        // prepare output
        let mut output = [E::Fr::zero(); RATE];
        for (o, s) in output.iter_mut().zip(state[..RATE].iter()) {
            *o = *s;
        }

        output
    }

    pub fn absorb_multiple(&mut self, input: &[E::Fr]) {
        // compute padding values
        let padding_strategy = DomainStrategy::CustomVariableLength;
        let padding_values = padding_strategy
            .generate_padding_values::<E>(input.len(), RATE);

        for inp in input.iter().chain(padding_values.iter()) {
            self.absorb(*inp)
        }
    }

    pub fn absorb(&mut self, input: E::Fr) {
        match self.mode {
            SpongeMode::Absorb(ref mut buf) => {
                // push value into buffer
                for el in buf.iter_mut() {
                    if el.is_none() {
                        // we still have empty room for values
                        *el = Some(input);
                        return;
                    }
                }

                // buffer is filled so unwrap them
                let mut unwrapped_buffer = [E::Fr::zero(); RATE];
                for (a, b) in unwrapped_buffer.iter_mut().zip(buf.iter_mut()) {
                    if let Some(val) = b {
                        *a = *val;
                        *b = None; // kind of resetting buffer
                    }
                }

                // here we can absorb values. run round function implicitly there
                absorb::<E, _, RATE, WIDTH>(&mut self.state, &mut unwrapped_buffer, self.params);

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
    }

    pub fn pad_if_necessary(&mut self) {
        match self.mode {
            SpongeMode::Absorb(ref mut buf) => {
                let unwrapped_buffer_len = buf.iter().filter(|el| el.is_some()).count();
                // compute padding values
                let padding_strategy = DomainStrategy::CustomVariableLength;
                let padding_values = padding_strategy
                    .generate_padding_values::<E>(unwrapped_buffer_len, RATE);
                let mut padding_values_it =  padding_values.iter().cloned();

                for b in buf {
                    if b.is_none() {
                        *b = padding_values_it.next()
                    }
                }
                assert!(padding_values_it.next().is_none());
            }
            SpongeMode::Squeeze(_) => (),
        }
    }

    pub fn squeeze(&mut self) -> Option<E::Fr> {
        loop {
            match self.mode {
                SpongeMode::Absorb(ref mut buf) => {
                    // buffer may not be filled fully so we may need padding.
                    let mut unwrapped_buffer = vec![];
                    for el in buf {
                        if let Some(value) = el {
                            unwrapped_buffer.push(*value);
                        } else {
                            // processing buffer was done and we need padding
                            break;
                        }
                    }
                    if unwrapped_buffer.is_empty() {
                        return None;
                    }

                    // make input array
                    let mut all_inputs = [E::Fr::zero(); RATE];
                    for (a, b) in all_inputs.iter_mut().zip(unwrapped_buffer) {
                        *a = b;
                    }

                    // permute state
                    absorb(&mut self.state, &all_inputs, self.params);

                    // push values into squeezing buffer for later squeezing
                    let mut squeeze_buffer = [None; RATE];
                    for (s, b) in self.state[..RATE].iter().zip(squeeze_buffer.iter_mut()) {
                        *b = Some(*s)
                    }

                    // we are switching squeezing mode so we can ignore to reset absorbing buffer
                    self.mode = SpongeMode::Squeeze(squeeze_buffer);
                }
                SpongeMode::Squeeze(ref mut buf) => {
                    for el in buf {
                        if let Some(value) = el.take() {
                            return Some(value);
                        }
                    }
                    return None;
                }
            };
        }
    }
}

fn absorb<E: Engine, P: HashParams<E, RATE, WIDTH>, const RATE: usize, const WIDTH: usize>(
    state: &mut [E::Fr; WIDTH],
    input: &[E::Fr; RATE],
    params: &P,
) {
    for (i, s) in input.iter().zip(state.iter_mut()) {
        s.add_assign(i);
    }
    generic_round_function(params, state, None);
}

pub fn generic_round_function<
    E: Engine,
    P: HashParams<E, RATE, WIDTH>,
    const RATE: usize,
    const WIDTH: usize,
>(
    params: &P,
    state: &mut [E::Fr; WIDTH],
    input: Option<[E::Fr; RATE]>,
) {
    if input.is_some() {
        unimplemented!("round function with absorb has not implemented yet");
    }

    match params.hash_family() {
        crate::traits::HashFamily::Rescue => {
            crate::rescue::rescue_round_function(params, state, input)
        }
        crate::traits::HashFamily::Poseidon => {
            crate::poseidon::poseidon_round_function(params, state, input)
        }
        crate::traits::HashFamily::RescuePrime => {
            crate::rescue_prime::rescue_prime_round_function(params, state, input)
        }
    }
}
