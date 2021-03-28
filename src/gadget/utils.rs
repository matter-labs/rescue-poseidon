use franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use franklin_crypto::plonk::circuit::linear_combination::LinearCombination;
use franklin_crypto::{
    bellman::{Engine, SynthesisError},
    plonk::circuit::allocated_num::Num,
};
use std::convert::TryInto;
// Multiplies matrix with a vector  and assigns result into same vector.
pub(crate) fn matrix_vector_product<E: Engine, CS: ConstraintSystem<E>, const DIM: usize>(
    cs: &mut CS,
    matrix: &[[E::Fr; DIM]; DIM],
    vector: &[LinearCombination<E>; DIM],
) -> Result<[LinearCombination<E>; DIM], SynthesisError> {
    let mut result: [LinearCombination<E>; DIM] = (0..DIM)
        .map(|_| LinearCombination::zero())
        .collect::<Vec<LinearCombination<E>>>()
        .try_into()
        .expect("vector of lc");
    // let mut result = [LinearCombination::zero(); DIM];
    let vec_as_nums = vector
        .iter()
        .map(|v| v.to_owned().into_num(cs).expect("into allocated num"))
        .collect::<Vec<Num<E>>>();
    for (i, matrix_row) in matrix.iter().enumerate() {
        // [fr, fr, fr] * [lc, lc, lc]
        for (coeff, num) in matrix_row.iter().zip(&vec_as_nums) {
            result[i].add_assign_number_with_coeff(&num, *coeff)
        }
    }

    Ok(result)
}

// Multiply sparse matrix and vector by exploiting sparsity of optimized matrixes.
pub(crate) fn mul_by_sparse_matrix<E: Engine, CS: ConstraintSystem<E>, const DIM: usize>(
    _cs: &mut CS,
    vector: &[LinearCombination<E>; DIM],
    matrix: &[[E::Fr; DIM]; DIM],
) -> [LinearCombination<E>; DIM] {
    assert_eq!(DIM, 3, "valid only for 3x3 matrix");
    let mut result: [LinearCombination<E>; DIM] = (0..DIM)
        .map(|_| LinearCombination::zero())
        .collect::<Vec<LinearCombination<E>>>()
        .try_into()
        .expect("vector of lc");

    for (a, b) in vector.iter().zip(matrix[0].iter()) {
        result[0].add_assign_scaled(a, *b);
    }

    result[1].add_assign_scaled(&vector[0], matrix[1][0]);
    result[1].add_assign(&vector[1]);

    result[2].add_assign_scaled(&vector[0], matrix[2][0]);
    result[2].add_assign(&vector[2]);

    result
}

#[cfg(test)]
mod test {
    use crate::tests::{init_cs, init_rng};
    use franklin_crypto::bellman::Field;
    use franklin_crypto::{
        bellman::pairing::bn256::{Bn256, Fr},
        plonk::circuit::{allocated_num::AllocatedNum, linear_combination::LinearCombination},
    };
    use rand::Rand;
    use std::convert::TryInto;

    #[test]
    fn test_matrix_product() {
        let cs = &mut init_cs::<Bn256>();
        let rng = &mut init_rng();

        const DIM: usize = 3;

        let mut vector_fe: [Fr; DIM] = [Fr::rand(rng); DIM];

        let mut vector_lc: [LinearCombination<_>; DIM] = (0..DIM)
            .map(|_| LinearCombination::zero())
            .collect::<Vec<LinearCombination<_>>>()
            .try_into()
            .expect("vector of lc");
        vector_fe
            .iter()
            .zip(vector_lc.iter_mut())
            .for_each(|(src, dst)| {
                *dst = LinearCombination::from(AllocatedNum::alloc(cs, || Ok(*src)).unwrap())
            });

        let mut matrix = [[Fr::zero(); DIM]; DIM];
        (0..9)
            .map(|_| Fr::rand(rng))
            .collect::<Vec<Fr>>()
            .chunks_exact(3)
            .zip(matrix.iter_mut())
            .for_each(|(src, dst)| *dst = src.try_into().expect("static vector"));

        matrix[1][1] = Fr::one();
        matrix[1][2] = Fr::zero();
        matrix[2][1] = Fr::zero();
        matrix[2][2] = Fr::one();

        crate::common::matrix::mmul_assign::<Bn256, DIM>(&matrix, &mut vector_fe);
        let actual = super::mul_by_sparse_matrix(cs, &vector_lc, &matrix);

        vector_fe.iter().zip(actual.iter()).for_each(|(fe, lc)| {
            let actual = lc.clone().into_num(cs).unwrap().get_value().unwrap();
            assert_eq!(*fe, actual);
        });
    }
}
