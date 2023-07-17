pub mod params;
pub mod poseidon2;

use self::params::Poseidon2Params;
pub use self::poseidon2::*;

use boojum::field::SmallField;
use boojum::cs::oracle::TreeHasher;
use franklin_crypto::bellman::{Engine, Field, PrimeField, PrimeFieldRepr};
use boojum::algebraic_props::round_function::AbsorptionModeTrait;

#[derive(Clone, Debug)]
pub struct Poseidon2<
    E: Engine,
    F: SmallField,
    M: AbsorptionModeTrait<E::Fr>,
    const RATE: usize,
    const WIDTH: usize
>{
    params: Poseidon2Params<E, RATE, WIDTH>,
    state: [E::Fr; WIDTH],
    buffer: [E::Fr; RATE],
    filled: usize,
    _marker: std::marker::PhantomData<(F, M)>,
}

impl<
    E: Engine,
    F: SmallField,
    M: AbsorptionModeTrait<E::Fr>,
    const RATE: usize,
    const WIDTH: usize,
> Poseidon2<E, F, M, RATE, WIDTH> {
    pub fn new() -> Self {
        assert!(Self::capasity_per_element() > 0);

        let params = Poseidon2Params::<E, RATE, WIDTH>::default();

        Self {
            params,
            state: [E::Fr::zero(); WIDTH],
            buffer: [E::Fr::zero(); RATE],
            filled: 0,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn capasity_per_element() -> usize {
        (E::Fr::CAPACITY as usize) / (F::CHAR_BITS as usize)
    }

    pub fn absorb_single(&mut self, value: &F) {
        let capasity_per_element = Self::capasity_per_element();
        debug_assert!(self.filled < RATE * capasity_per_element);
        let pos = self.filled / capasity_per_element;
        let exp = self.filled % capasity_per_element;

        let mut value_repr = <E::Fr as PrimeField>::Repr::from(value.as_u64());
        value_repr.shl((exp * capasity_per_element) as u32);

        self.buffer[pos].add_assign(&E::Fr::from_repr(value_repr).unwrap());
        self.filled += 1;

        if self.filled == RATE * capasity_per_element {
            for (dst, src) in self.state.iter_mut()
                .zip(self.buffer.iter_mut())
            {
                M::absorb(dst, src);
                *src = E::Fr::zero();
            }

            poseidon2_round_function(&mut self.state, &self.params);
            self.filled = 0;
        }
    }

    pub fn finalize(self) -> E::Fr {
        let Self {
            buffer,
            filled,
            state,
            params,
            ..
        } = self;

        let mut state = state;

        if filled > 0 {
            for (dst, src) in state.iter_mut().zip(buffer.iter()) {
                M::absorb(dst, src);
            }

            // We use just zeroes padding, so we don't need to call any padding function

            poseidon2_round_function(&mut state, &params);
        }

        state[0]
    }

    pub fn finalize_reset(&mut self) -> E::Fr {
        // reset
        let mut state = std::mem::replace(&mut self.state, [E::Fr::zero(); WIDTH]);
        let filled = self.filled;
        self.filled = 0;

        // run round function if necessary
        if filled > 0 {
            for (dst, src) in state.iter_mut().zip(self.buffer.iter_mut()) {
                M::absorb(dst, src);
                *src = E::Fr::zero();
            }

            // We use just zeroes padding, so we don't need to call any padding function

            poseidon2_round_function(&mut state, &self.params);
        }

        state[0]
    }
}

impl<
    E: Engine,
    F: SmallField,
    M: AbsorptionModeTrait<E::Fr>,
    const RATE: usize,
    const WIDTH: usize,
> TreeHasher<F> for Poseidon2<E, F, M, RATE, WIDTH> {
    type Output = E::Fr;

    #[inline]
    fn new() -> Self {
        Self::new()
    }

    #[inline]
    fn placeholder_output() -> Self::Output {
        E::Fr::zero()
    }

    #[inline]
    fn accumulate_into_leaf(&mut self, value: &F) {
        self.absorb_single(value);
    }

    #[inline]
    fn finalize_into_leaf_hash_and_reset(&mut self) -> Self::Output {
        self.finalize_reset()
    }

    #[inline]
    fn hash_into_leaf<'a, S: IntoIterator<Item = &'a F>>(source: S) -> Self::Output
    where
        F: 'a {
            let mut hasher = Self::new();

            for el in source.into_iter() {
                hasher.absorb_single(el);
            }
            hasher.finalize()
    }

    #[inline]
    fn hash_into_leaf_owned<S: IntoIterator<Item = F>>(source: S) -> Self::Output {
        let mut hasher = Self::new();

        for el in source.into_iter() {
            hasher.absorb_single(&el);
        }
        hasher.finalize()
    }

    #[inline]
    fn hash_into_node(left: &Self::Output, right: &Self::Output, _depth: usize) -> Self::Output {
        let params = Poseidon2Params::<E, RATE, WIDTH>::default();
        let mut state = [E::Fr::zero(); WIDTH];
        M::absorb(&mut state[0], left);
        M::absorb(&mut state[1], right);

        poseidon2_round_function(&mut state, &params);

        state[0]
    }
}
