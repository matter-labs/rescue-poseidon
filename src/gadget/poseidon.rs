use std::convert::TryInto;

use super::{sbox::*, utils::mul_by_sparse_matrix};
use super::{
    sponge::{
        GadgetSpongeMode, GadgetSpongePermutation, GadgetSpongeState, SpongeModes,
        StatefulSpongeGadget,
    },
    utils::matrix_vector_product,
};
use crate::{common::domain_strategy::DomainStrategy, sponge_gadget_impl, common::params::HasherParams};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::ConstraintSystem, plonk::circuit::boolean::Boolean,
};

use franklin_crypto::bellman::Field;
use franklin_crypto::bellman::SynthesisError;
use franklin_crypto::{
    bellman::Engine,
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter. Uses state-width=3 and rate=2.
pub fn poseidon_gadget<E, CS>(cs: &mut CS, input: &[Num<E>]) -> Result<[Num<E>; 2], SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    inner_poseidon_gadget::<_, _>(cs, input, DomainStrategy::CustomFixedLength)
}

/// Receives inputs whose length `unknown` prior (variable-length).
/// Also uses custom domain strategy which does not touch to value of capacity element
/// and does not apply any padding rule. Uses state-width=3 and rate=2.
pub fn poseidon_gadget_var_length<E, CS>(
    cs: &mut CS,
    input: &[Num<E>],
) -> Result<[Num<E>; 2], SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    inner_poseidon_gadget::<_, _>(cs, input, DomainStrategy::CustomVariableLength)
}

fn inner_poseidon_gadget<E, CS>(
    cs: &mut CS,
    input: &[Num<E>],
    domain: DomainStrategy<2>,
) -> Result<[Num<E>; 2], SynthesisError>
where
    E: Engine,
    CS: ConstraintSystem<E>,
{
    const STATE_WIDTH: usize = 3;
    const RATE: usize = 2;

    let result =
        super::hash::generic_hash::<E, _, PoseidonGadget<E, STATE_WIDTH, RATE>, STATE_WIDTH, RATE>(
            cs, input, domain,
        )?;

    Ok(result.try_into().expect("fixed length array"))
}

pub struct PoseidonGadget<E: Engine, const S: usize, const R: usize> {
    state: [LinearCombination<E>; S],
    params: HasherParams<E, S, R>,
    optimized_round_constants: Vec<[E::Fr; S]>,
    optimized_mds_matrixes: ([[E::Fr; S]; S], Vec<[[E::Fr; S]; S]>),
    sponge_mode: SpongeModes,
}

impl<E: Engine, const S: usize, const R: usize> Default for PoseidonGadget<E, S, R> {
    fn default() -> Self {
        let (params, _, optimized_round_constants, optimized_mds_matrixes) =
            crate::poseidon::params::poseidon_light_params();
        let initial_state: [LinearCombination<E>; S] = (0..S)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<E>>>()
            .try_into()
            .expect("vector of lc");
        Self {
            state: initial_state,
            params,
            optimized_round_constants,
            optimized_mds_matrixes,
            sponge_mode: SpongeModes::Standard(false),
        }
    }
}

sponge_gadget_impl!(PoseidonGadget<E, S, R>);

impl<E: Engine, const S: usize, const R: usize> GadgetSpongePermutation<E>
    for PoseidonGadget<E, S, R>
{
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
            .to_vec();
        constants_for_partial_rounds.push([E::Fr::zero(); S]);
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
            // reduce gate cost: LC -> Num -> LC
            for state in self.state.iter_mut() {
                let num = state.clone().into_num(cs).expect("a num");
                *state = LinearCombination::from(num.get_variable());
            }
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
    use crate::tests::init_cs;
    #[test]
    fn test_poseidon_light_sponge_with_custom_gate() {
        const STATE_WIDTH: usize = 3;
        const RATE: usize = 2;
        // TODO
        let cs = &mut init_cs();

        let mut el = Fr::one();
        el.double();

        let input = vec![el; RATE];

        let input_as_num = input
            .iter()
            .map(|el| Num::Variable(AllocatedNum::alloc(cs, || Ok(*el)).unwrap()))
            .collect::<Vec<Num<Bn256>>>();

        let mut poseidon_light_gadget = PoseidonGadget::<_, STATE_WIDTH, RATE>::default();
        poseidon_light_gadget.absorb(cs, &input_as_num).unwrap();
        let gadget_output: Vec<Num<Bn256>> = poseidon_light_gadget.squeeze(cs, None).unwrap();
        cs.finalize();
        assert!(cs.is_satisfied());

        println!("number of gates {}", cs.n());
        println!("last step number {}", cs.get_current_step_number());

        // poseidon_light original
        let mut poseidon_light = crate::poseidon::PoseidonHasher::<Bn256, 3, 2>::default();
        poseidon_light.absorb(&input);
        let output = poseidon_light.squeeze(None);

        for (sponge, gadget) in output.iter().zip(gadget_output.iter()) {
            assert_eq!(gadget.get_value().unwrap(), *sponge);
        }
    }
}
