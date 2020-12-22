use franklin_crypto::bellman::pairing::ff::{Field, PrimeField};
use franklin_crypto::bellman::pairing::Engine;

// Substitution box is non-linear part of permutation function.
// It basically computes power of each element in the state.
// Usually value of alpha is either 5 or 3. We keep a generic 
// handler other values of alpha. 
pub(crate) fn sbox<E: Engine>(alpha: E::Fr, state: &mut [E::Fr]) {
    match alpha.into_repr().as_ref()[0] {
        5 => {
            for el in state.iter_mut() {
                let mut quad = *el;
                quad.square();
                quad.square();
                el.mul_assign(&quad);
            }
        },
        3 => {
            for el in state.iter_mut() {
                let mut quad = *el;
                quad.square();
                el.mul_assign(&quad);
            }
        },
         _ => {
            for el in state.iter_mut() {
                *el = el.pow(alpha.into_repr());
            }
         }
    }
}
