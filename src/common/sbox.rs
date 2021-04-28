use crate::traits::Sbox;
use franklin_crypto::bellman::pairing::ff::Field;
use franklin_crypto::bellman::pairing::Engine;

// Substitution box is non-linear part of permutation function.
// It basically computes power of each element in the state.
// Usually value of alpha is either 5 or 3. We keep a generic
// handler other values of alpha.
pub(crate) fn sbox<E: Engine>(power: &Sbox, state: &mut [E::Fr]) {
    match power {
        Sbox::Alpha(alpha) => sbox_alpha::<E>(alpha, state),
        Sbox::AlphaInverse(alpha_inv) => sbox_alpha_inv::<E>(alpha_inv, state),
    }
}

pub(crate) fn sbox_alpha<E: Engine>(alpha: &u64, state: &mut [E::Fr]) {
    match alpha {
        5 => {
            for el in state.iter_mut() {
                let mut quad = *el;
                quad.square();
                quad.square();
                el.mul_assign(&quad);
            }
        }
        3 => {
            for el in state.iter_mut() {
                let mut quad = *el;
                quad.square();
                el.mul_assign(&quad);
            }
        }
        _ => {
            for el in state.iter_mut() {
                *el = el.pow(&[*alpha]);
            }
        }
    }
}

pub(crate) fn sbox_alpha_inv<E: Engine>(alpha_inv: &[u64; 4], state: &mut [E::Fr]) {
    for el in state.iter_mut() {
        *el = el.pow(alpha_inv);
    }
}
