use franklin_crypto::{
    bellman::{plonk::better_better_cs::cs::ConstraintSystem, Engine, SynthesisError},
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};

pub trait SpongeGadget<E: Engine, const RATE: usize, const WIDTH: usize> {
    fn specialize(
        &mut self,
        capacity_value: Option<LinearCombination<E>>,
    ) -> Result<(), SynthesisError>;

    fn absorb<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        input: &[Num<E>],
    ) -> Result<(), SynthesisError>;

    fn squeeze<CS: ConstraintSystem<E>>(
        &mut self,
        cs: &mut CS,
        number_of_elems: Option<usize>,
    ) -> Result<Vec<Num<E>>, SynthesisError>;

    fn reset(&mut self);
}
