use crate::traits::{HashFamily, HashParams, Sponge};
use franklin_crypto::bellman::{Engine, Field};

#[derive(Clone, Debug)]
pub enum SpongeModes {
    // Standard mode is stateless
    Standard(bool),
    // Duplex is statefull and maximum number of element "l" one can request
    // is equal to rate parameter.
    Duplex(bool),
}

#[derive(Debug, Clone)]
pub struct GenericSponge<
    'a,
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    const STATE_WIDTH: usize,
    const RATE: usize,
> {
    params: &'a P,
    state: [E::Fr; STATE_WIDTH],
    mode: SpongeModes,
}

impl<
        'a,
        E: Engine,
        P: HashParams<E, STATE_WIDTH, RATE>,
        const STATE_WIDTH: usize,
        const RATE: usize,
    > GenericSponge<'a, E, P, STATE_WIDTH, RATE>
{
    pub fn from_params(params: &'a P) -> Self {
        Self {
            params,
            state: [E::Fr::zero(); STATE_WIDTH],
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
    > Sponge<E, STATE_WIDTH, RATE> for GenericSponge<'a, E, P, STATE_WIDTH, RATE>
{
    fn specialize(&mut self, capacity_value: Option<E::Fr>) {
        let value = capacity_value.unwrap_or(E::Fr::zero());
        if let Some(last_el) = self.state.last_mut() {
            *last_el = value
        }
    }

    fn absorb(&mut self, input: &[E::Fr]) {
        assert!(!input.is_empty());
        let rate = RATE;

        match self.mode {
            SpongeModes::Standard(ref mut is_absorbed) => {
                assert_eq!(
                    input.len() % rate,
                    0,
                    "input length is not multiple of rate"
                );
                assert!(!*is_absorbed, "Sponge should be in in absorbtion phase");
                for elems in input.chunks_exact(rate) {
                    for (el, state) in elems.iter().zip(self.state.iter_mut()) {
                        state.add_assign(el);
                    }
                    generic_round_function(self.params, &mut self.state);
                    *is_absorbed = true;
                }
            }
            SpongeModes::Duplex(ref mut is_absorbed) => {
                assert!(!*is_absorbed, "Sponge should be in in absorbtion phase");
                assert!(
                    input.len() <= rate,
                    "duplex sponge can absorb max rate elems"
                );
                for (el, state) in input.iter().zip(self.state.iter_mut()) {
                    state.add_assign(el);
                }
                generic_round_function(self.params, &mut self.state);
                *is_absorbed = true;
            }
        }
    }

    fn squeeze(&mut self, number_of_elems: Option<usize>) -> Vec<E::Fr> {
        let rate = RATE;

        let mut out = vec![];

        match self.mode {
            SpongeModes::Standard(ref mut is_absorbed) => {
                assert!(*is_absorbed, "Sponge should be in in squeezing phase");
                if let Some(number_of_elems) = number_of_elems {
                    if number_of_elems <= rate {
                        out.extend_from_slice(&self.state[..number_of_elems]);
                    } else {
                        let original_number_of_elems = number_of_elems;

                        let number_of_iters = if number_of_elems % rate != 0 {
                            (number_of_elems + (rate - (number_of_elems % rate))) / rate
                        } else {
                            number_of_elems / rate
                        };

                        for _ in 0..number_of_iters {
                            out.extend_from_slice(&self.state[..rate]);
                            generic_round_function(self.params, &mut self.state);
                        }

                        out.truncate(original_number_of_elems);
                    }
                } else {
                    out.extend_from_slice(&self.state[..rate]);
                }
                *is_absorbed = false;
                self.reset();
                out
            }

            SpongeModes::Duplex(ref mut is_absorbed) => {
                assert!(*is_absorbed, "Sponge should be in in squeezing phase");
                let number_of_elems = if let Some(number_of_elems) = number_of_elems {
                    assert!(
                        number_of_elems <= rate,
                        "duplex sponge squeeze only as much as rate parameter"
                    );
                    number_of_elems
                } else {
                    rate
                };

                out.extend_from_slice(&self.state[..number_of_elems]);
                *is_absorbed = false;
                out
            }
        }
    }

    fn reset(&mut self) {
        self.state.iter_mut().for_each(|s| *s = E::Fr::zero());
    }
}

pub fn generic_round_function<
    E: Engine,
    P: HashParams<E, STATE_WIDTH, RATE>,
    const STATE_WIDTH: usize,
    const RATE: usize,
>(
    params: &P,
    state: &mut [E::Fr; STATE_WIDTH],
) {
    match params.hash_family() {
        HashFamily::Rescue => crate::rescue::rescue_round_function(params, state),
        HashFamily::Poseidon => crate::poseidon::poseidon_round_function(params, state),
        HashFamily::RescuePrime => crate::rescue_prime::rescue_prime_round_function(params, state),
    }
}
