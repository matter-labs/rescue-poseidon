use franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use franklin_crypto::plonk::circuit::linear_combination::LinearCombination;
use franklin_crypto::{
    bellman::{Engine, SynthesisError},
    plonk::circuit::allocated_num::Num,
};

// Multiplies matrix with a vector  and assigns result into same vector.
pub(crate) fn matrix_vector_product<E: Engine, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    matrix: &[Vec<E::Fr>],
    vector: &[LinearCombination<E>],
) -> Result<Vec<LinearCombination<E>>, SynthesisError> {
    let mut result = vec![LinearCombination::zero(); vector.len()];
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
pub(crate) fn mul_by_sparse_matrix<E: Engine, CS: ConstraintSystem<E>>(
    _cs: &mut CS,
    vector: &[LinearCombination<E>],
    matrix: &[Vec<E::Fr>],
) -> Vec<LinearCombination<E>> {
    let mut result = vec![LinearCombination::zero(); vector.len()];

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

    #[test]
    fn test_matrix_product() {
        let cs = &mut init_cs::<Bn256>();
        let rng = &mut init_rng();

        let mut vector_fe: Vec<Fr> = (0..3).map(|_| Fr::rand(rng)).collect();

        let vector_lc = vector_fe
            .iter()
            .map(|el| LinearCombination::from(AllocatedNum::alloc(cs, || Ok(*el)).unwrap()))
            .collect::<Vec<LinearCombination<_>>>();

        let mut matrix = (0..9)
            .map(|_| Fr::rand(rng))
            .collect::<Vec<Fr>>()
            .chunks_exact(3)
            .map(|c| c.to_vec())
            .collect::<Vec<Vec<Fr>>>();

        matrix[1][1] = Fr::one();
        matrix[1][2] = Fr::zero();
        matrix[2][1] = Fr::zero();
        matrix[2][2] = Fr::one();

        crate::common::matrix::mmul_assign::<Bn256>(&matrix, &mut vector_fe);
        let actual = super::mul_by_sparse_matrix(cs, &vector_lc, &matrix);

        vector_fe.iter().zip(actual).for_each(|(fe, lc)| {
            let actual = lc.clone().into_num(cs).unwrap().get_value().unwrap();
            assert_eq!(*fe, actual);
        });
    }
}
