use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;
// use poseidon_hash::StatefulSponge;

pub(crate) mod state;
pub trait SpongeParams {
    fn rate(&self) -> usize;
}
pub trait SpongeState<E: Engine, const S: usize> {
    fn state_as_ref(&self) -> &[E::Fr; S];
    fn state_as_mut(&mut self) -> &mut [E::Fr; S];
}

pub trait SpongePermutation<E: Engine> {
    fn permutation(&mut self);
}
#[derive(Clone, Debug)]
pub enum SpongeModes {
    // Standard mode is stateless
    Standard(bool),
    // Duplex is statefull and maximum number of element "l" one can request
    // is equal to rate parameter.
    Duplex(bool),
}

pub trait SpongeMode<E: Engine> {
    fn get_mode(&self) -> SpongeModes;
    fn update_mode(&mut self, mode: SpongeModes);
    fn mode_as_mut(&mut self) -> &mut SpongeModes;
}

pub trait StatefulSponge<E: Engine, const S: usize, const R: usize>:
    SpongeState<E, S> + SpongePermutation<E> + SpongeMode<E> + Default
{
    fn specialize(&mut self, capacity_value: Option<E::Fr>) {
        let state = self.state_as_mut();
        let value = capacity_value.unwrap_or(E::Fr::zero());
        if let Some(last_el) = state.last_mut() {
            *last_el = value
        }
    }

    /// In the absorbing phase, the r-element input blocks are summed into the first r
    /// elements of the state, interleaved with applications of the  permutation.
    /// When all input blocks are processed, the sponge construction switches to the
    /// squeezing phase.
    fn absorb(&mut self, input: &[E::Fr]) {
        assert!(!input.is_empty());
        let rate = R;

        match self.get_mode() {
            // TODO
            SpongeModes::Standard(is_absorbed) => {
                assert_eq!(
                    input.len() % rate,
                    0,
                    "input length is not multiple of rate"
                );
                assert!(!is_absorbed, "Sponge should be in in absorbtion phase");
                for elems in input.chunks_exact(rate) {
                    for (el, state) in elems.iter().zip(self.state_as_mut().iter_mut()) {
                        state.add_assign(el);
                    }
                    self.permutation();
                    // *is_absorbed = true; // absorbed
                    self.update_mode(SpongeModes::Standard(true));
                }
            }
            SpongeModes::Duplex(is_absorbed) => {
                assert!(!is_absorbed, "Sponge should be in in absorbtion phase");
                assert!(
                    input.len() <= rate,
                    "duplex sponge can absorb max rate elems"
                );
                // If state already squeezed then discard buffer. We don't need to
                // accumulate any value here because we alread stored in top of function
                // TODO
                for (el, state) in input.iter().zip(self.state_as_mut().iter_mut()) {
                    state.add_assign(el);
                }
                self.permutation();
                self.update_mode(SpongeModes::Standard(true));
            }
        }
    }

    /// In the squeezing phase, the first r-elements of the state are
    /// returned as output blocks, interleaved with applications of the permutation.
    /// The number of output blocks is chosen at will by the user.
    fn squeeze(&mut self, number_of_elems: Option<usize>) -> Vec<E::Fr> {
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
                            self.permutation();
                        }

                        out.truncate(original_number_of_elems);
                    }
                } else {
                    out.extend_from_slice(&self.state_as_ref()[..rate]);
                }
                self.update_mode(SpongeModes::Standard(false));
                self.reset();
                out
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
                out
            }
        }
    }

    fn reset(&mut self) {
        self.state_as_mut()
            .iter_mut()
            .for_each(|s| *s = E::Fr::zero());
    }
}
