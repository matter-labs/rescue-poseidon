use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;

pub(crate) mod state;
pub trait SpongeParams {
    fn rate(&self) -> usize;
}
pub trait SpongeState<E: Engine> {
    fn state_as_ref(&self) -> &[E::Fr];
    fn storage_as_ref(&self) -> &[E::Fr];
    fn state_as_mut(&mut self) -> &mut [E::Fr];
    fn storage_as_mut(&mut self) -> &mut Vec<E::Fr>;
}

pub trait SpongePermutation<E: Engine> {
    fn permutation(&mut self);
}
#[derive(Clone, Debug)]
pub enum SpongeModes<E: Engine> {
    // Standard mode is stateless
    Standard,
    // Duplex is statefull and maximum number of element "l" one can request
    // is equal to rate parameter.
    Duplex(Vec<E::Fr>),
}

pub trait SpongeMode<E: Engine> {
    fn get_mode(&self) -> SpongeModes<E>;
    fn update_mode(&mut self, mode: SpongeModes<E>);
}

/// In the absorbing phase, the r-element input blocks are summed into the first r
/// elements of the state, interleaved with applications of the  permutation.
/// When all input blocks are processed, the sponge construction switches to the
/// squeezing phase. In the squeezing phase, the first r-elements of the state are
/// returned as output blocks, interleaved with applications of the permutation.
/// The number of output blocks is chosen at will by the user.
pub trait StatefulSponge<E: Engine>:
    SpongeState<E> + SpongePermutation<E> + SpongeParams + SpongeMode<E> + Default
{
    fn specialize(&mut self, capacity_value: Option<E::Fr>) {
        let state = self.state_as_mut();
        let value = capacity_value.unwrap_or(E::Fr::zero());
        if let Some(last_el) = state.last_mut() {
            *last_el = value
        }
    }

    fn absorb(&mut self, input: E::Fr) {
        let rate = self.rate();
        if self.storage_as_ref().len() < rate {
            self.storage_as_mut().push(input);
        }

        match self.get_mode() {
            SpongeModes::Standard => {
                if self.storage_as_ref().len() == rate {                    
                    permute(self);
                }
            }
            SpongeModes::Duplex(_) => {
                // If state already squeezed then discard buffer. We don't need to
                // accumulate any value here because we alread stored in top of function
                self.update_mode(SpongeModes::Duplex(Vec::with_capacity(self.rate())));
            }
        }
    }

    fn absorb_multi(&mut self, input: &[E::Fr]) {
        let rate = self.rate();
        assert!(!input.is_empty());
        assert_eq!(input.len() % rate, 0);
        for i in 0..input.len() / rate {
            for value in input.iter().skip(i * rate).take(rate) {
                self.absorb(*value);
            }
        }
    }

    fn squeeze(&mut self) -> Vec<E::Fr> {
        assert!(
            self.storage_as_ref().is_empty(),
            "storage should be empty which also means permutation happened"
        );
        assert!(
            self.state_as_ref().iter().all(|el| !el.is_zero()),
            "state elements should not equal to zero"
        );

        self.state_as_ref()[..self.rate()].to_vec()
    }

    fn squeeze_single(&mut self) -> E::Fr {
        match self.get_mode() {
            SpongeModes::Standard => {
                unimplemented!("Only duplex sponge can squeeze single element")
            }
            SpongeModes::Duplex(mut buffer) => {
                // use already squeezed values
                if buffer.len() > 0 {
                    let out = buffer.remove(0);
                    self.update_mode(SpongeModes::Duplex(buffer));
                    return out;
                }
                // at least one element should have been absorbed
                assert!(self.storage_as_ref().len() >= 1);
                // pad elements in order to run a permutation
                // TODO:  PaddingStrategy enum can be used here
                while self.storage_as_ref().len() % self.rate() != 0 {
                    self.storage_as_mut().push(E::Fr::one());
                }
                permute(self);
                // squeeze from state
                buffer.extend_from_slice(&self.state_as_ref()[..self.rate()]);
                // output single element
                let out = buffer.remove(0);
                self.update_mode(SpongeModes::Duplex(buffer));
                out
            }
        }
    }

    fn squeeze_multi(&mut self, len: usize) -> Vec<E::Fr> {
        assert!(self.storage_as_ref().is_empty(), "storage should be empty");
        assert!(
            self.state_as_ref().iter().all(|el| !el.is_zero()),
            "state elements should not equal to zero"
        );
        assert!(
            len > 2,
            "length of requested output should be greater than rate param"
        );
        let mut output = vec![];

        while output.len() < len {
            for i in 0..self.rate() {
                output.push(self.state_as_ref()[i]);
            }

            if output.len() < len {
                self.permutation();
                self.storage_as_mut().truncate(0);
            }
        }
        assert!(output.len() == len);

        output
    }

    fn reset(&mut self) {
        self.storage_as_mut().truncate(0);
        self.state_as_mut()
            .iter_mut()
            .for_each(|s| *s = E::Fr::zero());
    }
}

fn permute<E: Engine, S: StatefulSponge<E>>(sponge: &mut S) {
    let storage_values = sponge.storage_as_ref().to_vec();
    for (value, state) in storage_values.iter().zip(sponge.state_as_mut().iter_mut()) {
        state.add_assign(value)
    }

    sponge.permutation();
    sponge.storage_as_mut().truncate(0);
}
