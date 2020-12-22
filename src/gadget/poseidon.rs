use super::{sbox::*, utils::mul_by_sparse_matrix};
use super::{
    sponge::{
        GadgetSpongeParams, GadgetSpongePermutation, GadgetSpongeState, StatefulSpongeGadget, GadgetSpongeMode, SpongeModes
    },
    utils::matrix_vector_product,
};
use crate::{common::padding::PaddingStrategy, sponge_gadget_impl, HasherParams};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::boolean::Boolean,
};

use franklin_crypto::bellman::Field;
use franklin_crypto::bellman::SynthesisError;
use franklin_crypto::{
    bellman::Engine,
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};

/// Constant-Input-Length Hashing. 
/// The capacity value is length x (2^64 ) + (o − 1) where o the output length. 
/// The padding consists of the field elements being 0.
pub fn poseidon_gadget_fixed_length<E, CS>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<Vec<Num<E>>, SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    super::hash::generic_hash::<E, _, PoseidonGadget<E>>(
        cs,
        input,
        PaddingStrategy::FixedLength(input.len()),
    )
}

/// Variable-Input-Length Hashing. 
/// The capacity value is 2^64 + (o − 1) where o the output length. 
/// The padding consists of one field element being 1, and the remaining elements being 0.
///  padding is necessary for variable-length inputs, even if the input is already (without delimiter) a multiple of the rate in length.
pub fn poseidon_gadget_var_length<E, CS>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<Vec<Num<E>>, SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    super::hash::generic_hash::<E, _, PoseidonGadget<E>>(
        cs,
        input,
        PaddingStrategy::VariableLength(input.len()),
    )
}

/// Similar to function with variable length input but with a small difference.
/// This function uses custom specialization with custom padding strategy.
pub fn poseidon_gadget<E, CS>(cs: &mut CS, input: &[Num<E>]) -> Result<Vec<Num<E>>, SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    super::hash::generic_hash::<E, _, PoseidonGadget<E>>(
        cs,
        input,
        PaddingStrategy::Custom(input.len()),
    )
}
/// Stateful poseidon 
pub struct PoseidonGadget<E: Engine> {
    state: Vec<LinearCombination<E>>,
    params: HasherParams<E>,
    optimized_round_constants: Vec<Vec<E::Fr>>,
    optimized_mds_matrixes: (Vec<Vec<E::Fr>>, Vec<Vec<Vec<E::Fr>>>),
    sponge_mode: SpongeModes<E>,
    tmp_storage: Vec<Num<E>>,
}

impl<E: Engine> Default for PoseidonGadget<E> {
    fn default() -> Self {
        let (params, _, optimized_round_constants, optimized_mds_matrixes) =
            crate::poseidon::params::poseidon_light_params();
        Self {
            state: vec![LinearCombination::zero(); params.state_width],
            tmp_storage: Vec::with_capacity(params.rate),
            params,
            optimized_round_constants,
            optimized_mds_matrixes,
            sponge_mode: SpongeModes::Standard,
        }
    }
}

sponge_gadget_impl!(PoseidonGadget<E>);

// permutation happens in 4 full, 33 partial and 4 full rounds consecutively
// total cost 2 + 3*2 + 8*3*(2+2) = 104
impl<E: Engine> GadgetSpongePermutation<E> for PoseidonGadget<E> {
    fn permutation<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        _should_permute: &Boolean,
    ) -> Result<(), SynthesisError> {
        debug_assert!(self.params.full_rounds % 2 == 0);

        let half_of_full_rounds = self.params.full_rounds / 2;

        let (m_prime, sparse_matrixes) = &self.optimized_mds_matrixes;
        let optimized_round_constants = &self.optimized_round_constants;

        // first full rounds
        for round in 0..half_of_full_rounds {
            let round_constants = &optimized_round_constants[round];

            // add round constatnts
            for (s, c) in self.state.iter_mut().zip(round_constants.iter()) {
                s.add_assign_constant(*c);
            }
            // non linear sbox
            sbox_quintic::<E, _>(cs, &mut self.state)?;

            // mul state by mds
            self.state = matrix_vector_product(cs, &self.params.mds_matrix, &self.state)?;
        }

        self.state
            .iter_mut()
            .zip(optimized_round_constants[half_of_full_rounds].iter())
            .for_each(|(a, b)| a.add_assign_constant(*b));

        self.state = matrix_vector_product(cs, &m_prime, &self.state)?;

        let mut constants_for_partial_rounds = optimized_round_constants
            [half_of_full_rounds + 1..half_of_full_rounds + self.params.partial_rounds]
            // TODOC
            .to_vec();
        constants_for_partial_rounds.push(vec![E::Fr::zero(); 3]);
        // in order to reduce gate number we merge two consecutive iteration
        // which costs 2 gates per each
        for (round_constant, sparse_matrix) in constants_for_partial_rounds
            [..constants_for_partial_rounds.len() - 1]
            .chunks(2)
            .zip(sparse_matrixes[..sparse_matrixes.len() - 1].chunks(2))
        {
            // first
            sbox_quintic::<E, _>(cs, &mut self.state[..1])?;
            self.state[0].add_assign_constant(round_constant[0][0]);
            self.state = mul_by_sparse_matrix(cs, &self.state, &sparse_matrix[0]);

            // second
            sbox_quintic::<E, _>(cs, &mut self.state[..1])?;
            self.state[0].add_assign_constant(round_constant[1][0]);
            self.state = mul_by_sparse_matrix(cs, &self.state, &sparse_matrix[1]);
            self.state = lc_to_num_to_lc(cs, &self.state)?;
        }

        sbox_quintic::<E, _>(cs, &mut self.state[..1])?;
        self.state[0].add_assign_constant(constants_for_partial_rounds.last().unwrap()[0]);
        self.state = mul_by_sparse_matrix(cs, &self.state, &sparse_matrixes.last().unwrap());

        // second full round
        for round in (self.params.partial_rounds + half_of_full_rounds)
            ..(self.params.partial_rounds + self.params.full_rounds)
        {
            let round_constants = &optimized_round_constants[round];

            // add round constatnts
            for (s, c) in self.state.iter_mut().zip(round_constants.iter()) {
                s.add_assign_constant(*c);
            }

            sbox_quintic::<E, _>(cs, &mut self.state)?;

            // mul state by mds
            self.state = matrix_vector_product(cs, &self.params.mds_matrix, &self.state)?;
        }

        Ok(())
    }
}

fn lc_to_num_to_lc<E: Engine, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    input: &[LinearCombination<E>],
) -> Result<Vec<LinearCombination<E>>, SynthesisError> {
    Ok(input
        .iter()
        .cloned()
        .map(|lc| lc.into_num(cs).unwrap())
        .map(|num| LinearCombination::from(num.get_variable()))
        .collect::<Vec<LinearCombination<E>>>())
}

#[cfg(test)]
mod test {
    use franklin_crypto::bellman::{
        pairing::bn256::{Bn256, Fr},
        plonk::better_better_cs::cs::ConstraintSystem,
    };
    use franklin_crypto::{
        bellman::Field,
        plonk::circuit::allocated_num::{AllocatedNum, Num},
    };

    use crate::gadget::sponge::StatefulSpongeGadget;
    use crate::sponge::StatefulSponge;

    use super::PoseidonGadget;
    use crate::tests::{init_cs};
    #[test]
    fn test_poseidon_light_sponge_with_custom_gate() {
        let cs = &mut init_cs();

        let mut el = Fr::one();
        el.double();

        let input = vec![el; 2];

        let input_as_num = input
            .iter()
            .map(|el| Num::Variable(AllocatedNum::alloc(cs, || Ok(*el)).unwrap()))
            .collect::<Vec<Num<Bn256>>>();

        let mut poseidon_light_gadget = PoseidonGadget::default();
        poseidon_light_gadget
            .absorb_multi(cs, &input_as_num)
            .unwrap();
        let gadget_output = poseidon_light_gadget.squeeze(cs).unwrap();
        cs.finalize();
        assert!(cs.is_satisfied());

        println!("number of gates {}", cs.n());
        println!("last step number {}", cs.get_current_step_number());

        // poseidon_light original
        let mut poseidon_light = crate::poseidon::PoseidonHasher::<Bn256>::default();
        // let mut poseidon_light = crate::poseidon::PoseidonHasher::<Bn256>::default();
        poseidon_light.absorb_multi(&input);
        let output = poseidon_light.squeeze();

        for (sponge, gadget) in output.iter().zip(gadget_output.iter()) {
            assert_eq!(gadget.get_value().unwrap(), *sponge);
        }
    }
}
