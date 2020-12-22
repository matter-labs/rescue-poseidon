use franklin_crypto::bellman::{Engine, Field};
use std::ops::Range;

// We can reduce cost of each partial round by using an optimization from
// original paper. Appendix-B explains details.
pub(crate) fn compute_optimized_matrixes<E: Engine>(
    number_of_rounds: usize,
    original_mds: &[Vec<E::Fr>],
) -> (Vec<Vec<E::Fr>>, Vec<Vec<Vec<E::Fr>>>) {
    let state_width = original_mds.len();
    let original_mds = transpose::<E>(original_mds);
    let mut matrix = original_mds.to_vec();
    let mut m_prime = identity::<E>(state_width);
    let mut sparse_matrixes = vec![];
    for _ in 0..number_of_rounds {
        // M'
        let m_hat = sub_matrix::<E>(&matrix, 1..state_width, 1..state_width);
        m_prime = identity::<E>(state_width);
        set_sub_matrix::<E>(&mut m_prime, 1..state_width, 1..state_width, &m_hat);

        // M"
        let w = sub_matrix::<E>(&matrix, 1..state_width, 0..1);
        let v = sub_matrix::<E>(&matrix, 0..1, 1..state_width);

        let m_hat_inv = try_inverse::<E>(&m_hat).expect("inverse");
        let w_hat = multiply::<E>(&m_hat_inv, &w);

        let mut sparse_matrix = identity::<E>(state_width);
        sparse_matrix[0][0] = matrix[0][0];
        set_sub_matrix::<E>(&mut sparse_matrix, 0..1, 1..state_width, &v);
        set_sub_matrix::<E>(&mut sparse_matrix, 1..state_width, 0..1, &w_hat);
        {
            // sanity check
            let actual = multiply::<E>(&m_prime, &sparse_matrix);
            assert_eq!(matrix, actual);
        }
        sparse_matrixes.push(transpose::<E>(&sparse_matrix));
        matrix = multiply::<E>(&original_mds, &m_prime);
    }

    sparse_matrixes.reverse();
    sparse_matrixes.iter().chain(&[m_prime.clone()]).for_each(|matrix| {
        let _ = try_inverse::<E>(matrix).expect("should have matrixz");
    });

    (transpose::<E>(&m_prime), sparse_matrixes)
}

// Multiply sparse matrix and vector by exploiting sparsity of optimized matrixes.
pub(crate) fn mul_by_sparse_matrix<E: Engine>(
    vector: &[E::Fr],
    matrix: &[Vec<E::Fr>], // transposed
) -> Vec<E::Fr> {
    let mut result = vec![E::Fr::zero(); vector.len()];

    for (a, b) in vector.iter().zip(matrix[0].iter()) {
        let mut tmp = a.clone();
        tmp.mul_assign(&b);
        result[0].add_assign(&tmp);
    }

    let mut tmp = matrix[1][0];
    tmp.mul_assign(&vector[0]);
    tmp.add_assign(&vector[1]);
    result[1] = tmp;

    let mut tmp = matrix[2][0];
    tmp.mul_assign(&vector[0]);
    tmp.add_assign(&vector[2]);
    result[2] = tmp;

    result
}

// Decontructs a sub matrix
pub(crate) fn sub_matrix<E: Engine>(
    matrix: &[Vec<E::Fr>],
    row_range: std::ops::Range<usize>,
    col_range: std::ops::Range<usize>,
) -> Vec<Vec<E::Fr>> {
    let mut values = matrix.to_vec();
    let values: Vec<Vec<E::Fr>> = values
        .drain(row_range)
        .collect::<Vec<Vec<E::Fr>>>()
        .iter_mut()
        .map(|row| row.drain(col_range.clone()).collect())
        .collect();

    values
}

// Injects a lower dimension matrix into higher one.
pub(crate) fn set_sub_matrix<E: Engine>(
    matrix: &mut [Vec<E::Fr>],
    row_range: Range<usize>,
    col_range: Range<usize>,
    sub_matrix: &[Vec<E::Fr>],
) {
    for (row_a, row_b) in matrix[row_range].iter_mut().zip(sub_matrix.iter()) {
        for (col_a, col_b) in row_a[col_range.clone()].iter_mut().zip(row_b.iter()) {
            *col_a = *col_b;
        }
    }
}

// Multiplies matrix with a vector  and assigns result into same vector.
pub(crate) fn mmul_assign<E: Engine>(matrix: &[Vec<E::Fr>], vector: &mut [E::Fr]) {
    // [M]xv
    assert!(!matrix.is_empty());
    assert!(!vector.is_empty());
    let row_len = matrix.len();
    assert_eq!(row_len, vector.len());

    let mut result = vec![E::Fr::zero(); row_len];
    let row_len = matrix[0].len();
    matrix.iter().for_each(|row| assert_eq!(row_len, row.len()));
    for col in 0..row_len {
        result[col] = crate::common::utils::scalar_product::<E>(vector, &matrix[col]);
    }
    vector.copy_from_slice(&result[..]);
}

// Multiplies two same dimension matrixes.
pub(crate) fn multiply<E: Engine>(m1: &[Vec<E::Fr>], m2: &[Vec<E::Fr>]) -> Vec<Vec<E::Fr>> {
    assert_eq!(m1[0].len(), m2.len());
    let number_of_rows = m1.len();
    let number_of_cols = m2[0].len();
    let transposed_m2 = transpose::<E>(m2);

    let mut result = vec![vec![E::Fr::zero(); number_of_cols]; number_of_rows];

    for (i, rv) in m1.iter().enumerate() {
        for (j, cv) in transposed_m2.iter().enumerate() {
            result[i][j] = crate::common::utils::scalar_product::<E>(rv, cv);
        }
    }

    result
}
// Transpose of a matrix.
pub(crate) fn transpose<E: Engine>(matrix: &[Vec<E::Fr>]) -> Vec<Vec<E::Fr>> {
    let row_len = matrix.len();
    let col_len = matrix[0].len();

    let mut values = vec![vec![E::Fr::zero(); row_len]; col_len];
    for i in 0..row_len {
        for j in 0..col_len {
            values[j][i] = matrix[i][j];
        }
    }

    values
}

// Computes inverse of 2-d or 3-d matrixes.
pub(crate) fn try_inverse<E: Engine>(m: &[Vec<E::Fr>]) -> Option<Vec<Vec<E::Fr>>> {
    assert_eq!(m[0].len(), m.len());

    let dim = m.len();

    match dim {
        2 => try_inverse_dim_2::<E>(m),
        3 => try_inverse_dim_3::<E>(m),
        _ => unimplemented!(),
    }
}

// Computes inverse of 2x2 matrix.
fn try_inverse_dim_2<E: Engine>(m: &[Vec<E::Fr>]) -> Option<Vec<Vec<E::Fr>>> {
    let determinant = {
        let mut a = m[0][0];
        a.mul_assign(&m[1][1]);

        let mut b = m[1][0];
        b.mul_assign(&m[0][1]);

        a.sub_assign(&b);

        a
    };

    // if determinant.is_zero() {
    //     return None;
    // }
    let mut result = vec![vec![E::Fr::zero(); 2]; 2];
    let det_inv = if let Some(inv) = determinant.inverse() {
        inv
    } else {
        return None;
    };

    // m22 / determinant;
    result[0][0] = {
        let mut tmp = m[1][1];
        tmp.mul_assign(&det_inv);
        tmp
    };
    // -m12 / determinant;
    result[0][1] = {
        let mut tmp = m[0][1];
        tmp.negate();
        tmp.mul_assign(&det_inv);
        tmp
    };
    // -m21 / determinant;
    result[1][0] = {
        let mut tmp = m[1][0];
        tmp.negate();
        tmp.mul_assign(&det_inv);
        tmp
    };
    // m11 / determinant;
    result[1][1] = {
        let mut tmp = m[0][0];
        tmp.mul_assign(&det_inv);
        tmp
    };

    Some(result)
}

// Computes inverse of 3x3 matrix.
fn try_inverse_dim_3<E: Engine>(m: &[Vec<E::Fr>]) -> Option<Vec<Vec<E::Fr>>> {
    // m22 * m33 - m32 * m23;
    let minor_m12_m23 = {
        let mut a = m[1][1];
        a.mul_assign(&m[2][2]);

        let mut b = m[2][1];
        b.mul_assign(&m[1][2]);

        a.sub_assign(&b);

        a
    };

    //  m21 * m33 - m31 * m23;
    let minor_m11_m23 = {
        let mut a = m[1][0];
        a.mul_assign(&m[2][2]);

        let mut b = m[2][0];
        b.mul_assign(&m[1][2]);

        a.sub_assign(&b);

        a
    };
    // m21 * m32 - m31 * m22;
    let minor_m11_m22 = {
        let mut a = m[1][0];
        a.mul_assign(&m[2][1]);

        let mut b = m[2][0];
        b.mul_assign(&m[1][1]);

        a.sub_assign(&b);

        a
    };

    // m11 * minor_m12_m23 - m12 * minor_m11_m23 + m13 * minor_m11_m22;
    let determinant = {
        let mut a = m[0][0];
        a.mul_assign(&minor_m12_m23);

        let mut b = m[0][1];
        b.mul_assign(&minor_m11_m23);

        let mut c = m[0][2];
        c.mul_assign(&minor_m11_m22);

        a.sub_assign(&b);
        a.add_assign(&c);

        a
    };

    if determinant.is_zero() {
        // matrix is not invertible
        return None;
    }

    let mut result = vec![vec![E::Fr::zero(); m[0].len()]; m.len()];
    let det_inv = if let Some(inv) = determinant.inverse() {
        inv
    } else {
        return None;
    };
    //  minor_m12_m23 / determinant;
    result[0][0] = {
        let mut tmp = minor_m12_m23.clone();
        tmp.mul_assign(&det_inv);

        tmp
    };

    // (m13 * m32 - m33 * m12) / determinant;
    result[0][1] = {
        let mut a = m[0][2];
        a.mul_assign(&m[2][1]);

        let mut b = m[2][2];
        b.mul_assign(&m[0][1]);

        a.sub_assign(&b);

        a.mul_assign(&det_inv);

        a
    };

    // (m12 * m23 - m22 * m13) / determinant;
    result[0][2] = {
        let mut a = m[0][1];
        a.mul_assign(&m[1][2]);

        let mut b = m[1][1];
        b.mul_assign(&m[0][2]);

        a.sub_assign(&b);

        a.mul_assign(&det_inv);

        a
    };
    // -minor_m11_m23 / determinant;
    result[1][0] = {
        let mut tmp = minor_m11_m23;
        tmp.negate();
        tmp.mul_assign(&det_inv);

        tmp
    };
    // (m11 * m33 - m31 * m13) / determinant;
    result[1][1] = {
        let mut a = m[0][0];
        a.mul_assign(&m[2][2]);

        let mut b = m[2][0];
        b.mul_assign(&m[0][2]);

        a.sub_assign(&b);

        a.mul_assign(&det_inv);

        a
    };
    // (m13 * m21 - m23 * m11) / determinant;
    result[1][2] = {
        let mut a = m[0][2];
        a.mul_assign(&m[1][0]);

        let mut b = m[1][2];
        b.mul_assign(&m[0][0]);

        a.sub_assign(&b);

        a.mul_assign(&det_inv);

        a
    };
    // minor_m11_m22 / determinant;
    result[2][0] = {
        let mut tmp = minor_m11_m22;
        tmp.mul_assign(&det_inv);

        tmp
    };
    // m12 * m31 - m32 * m11) / determinant;
    result[2][1] = {
        let mut a = m[0][1];
        a.mul_assign(&m[2][0]);

        let mut b = m[2][1];
        b.mul_assign(&m[0][0]);

        a.sub_assign(&b);

        a.mul_assign(&det_inv);

        a
    };
    // (m11 * m22 - m21 * m12) / determinant;
    result[2][2] = {
        let mut a = m[0][0];
        a.mul_assign(&m[1][1]);

        let mut b = m[1][0];
        b.mul_assign(&m[0][1]);

        a.sub_assign(&b);

        a.mul_assign(&det_inv);

        a
    };

    Some(result)
}

// Computes identity of given dimension.
fn identity<E: Engine>(dimension: usize) -> Vec<Vec<E::Fr>> {
    let mut identity = vec![vec![E::Fr::zero(); dimension]; dimension];
    for i in 0..dimension {
        for j in 0..dimension {
            let el = if i == j { E::Fr::one() } else { E::Fr::zero() };
            identity[i][j] = el;
        }
    }

    identity
}

#[cfg(test)]
mod test {
    use crate::tests::init_rng;

    use super::*;
    use franklin_crypto::bellman::bn256::{Bn256, Fr};
    use franklin_crypto::bellman::PrimeField;
    use rand::Rand;
    #[test]
    fn test_matrix_inverese() {
        let one = Fr::one();
        let mut two = one.clone();
        two.add_assign(&one);
        let mut three = two.clone();
        three.add_assign(&one);

        let dim = 3;
        let values = vec![
            vec![two, one, one],
            vec![three, two, one],
            vec![two, one, two],
        ];

        let _ = try_inverse::<Bn256>(&values);

        assert_eq!(
            identity::<Bn256>(dim),
            multiply::<Bn256>(&try_inverse::<Bn256>(&values).expect("inverse"), &values)
        );
    }

    #[test]
    fn test_matrix_deconstruction() {
        let one = Fr::one();
        let mut two = one.clone();
        two.add_assign(&one);
        let mut three = two.clone();
        three.add_assign(&one);

        let matrix = vec![
            vec![two, one, one],
            vec![three, two, one],
            vec![two, one, two],
        ];

        {
            let expected = vec![vec![two, one], vec![one, two]];
            let actual = sub_matrix::<Bn256>(&matrix, 1..3, 1..3);
            assert_eq!(expected, actual);
        }
        {
            let expected = vec![vec![two], vec![three], vec![two]];
            let actual = sub_matrix::<Bn256>(&matrix, 0..3, 0..1);
            assert_eq!(expected, actual);
        }
    }

    #[test]
    fn test_inject_sub_matrix() {
        let zero = Fr::zero();
        let one = Fr::one();
        let mut two = one.clone();
        two.add_assign(&one);
        let mut three = two.clone();
        three.add_assign(&one);

        let mut matrix = vec![
            vec![two, one, one],
            vec![three, two, one],
            vec![two, one, two],
        ];
        let sub_matrix = vec![vec![zero, zero], vec![zero, zero]];

        let expected_matrix = vec![
            vec![two, one, one],
            vec![three, zero, zero],
            vec![two, zero, zero],
        ];

        set_sub_matrix::<Bn256>(&mut matrix, 1..3, 1..3, &sub_matrix);
        assert_eq!(expected_matrix, matrix);
    }

    #[test]
    fn test_matrix_transpose() {
        let rng = &mut init_rng();

        let dim = 3;
        for _ in 0..10 {
            let mut matrix = vec![vec![Fr::zero(); dim]; dim];
            for i in 0..dim {
                for j in 0..dim {
                    matrix[i][j] = Fr::rand(rng);
                }
            }
            assert_eq!(transpose::<Bn256>(&transpose::<Bn256>(&matrix)), matrix);
        }
    }

    #[test]
    fn test_optimized_matrixes() {
        use franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};

        let rng = &mut init_rng();

        let state_width = 3usize;

        let original_mds_1d =
            crate::common::utils::construct_mds_matrix::<Bn256, _>(state_width, rng);

        let original_mds = original_mds_1d
            .chunks_exact(state_width)
            .map(|els| els.to_vec())
            .collect::<Vec<Vec<Fr>>>();

        let (_, _) = compute_optimized_matrixes::<Bn256>(5, &original_mds);
    }

    #[test]
    fn test_sparse_multiplication() {
        // sparse matrix [[a, b, c], [d, 1, 0], [g, 0, 1]]
    }

    // #[bench]
    // fn bench_sparse_multiplication1(b: &mut Bencher) {
    //     let rng = &mut init_rng();
    //     // we need to bench sparse multiplication in two different ways
    //     // 1. use major row mode
    //     // 2. use single vector and chunks it by 3
    //     let dim = 3;
    //     let mut state = vec![Fr::zero(); dim];
    //     for i in 0..dim{
    //         state[i] = Fr::rand(rng);
    //     }

    //     // create random 3x3 matrix
    //     let mut two_d_matrix = vec![vec![Fr::zero(); dim]; dim];
    //     for i in 0..dim {
    //         for j in 0..dim {
    //             two_d_matrix[i][j] = Fr::rand(rng);
    //         }
    //     }

    //     b.iter(|| mul_by_sparse_matrix::<Bn256>(&state, &two_d_matrix));
    // }

    fn int_to_fe<E: Engine>(elements: &[i8]) -> Vec<E::Fr> {
        elements
            .iter()
            .map(|el| {
                let mut repr = <E::Fr as PrimeField>::Repr::default();
                repr.as_mut()[0] = *el as u64;
                let mut fe = E::Fr::from_repr(repr).expect("valid fe");
                if *el < 0 {
                    fe.negate();
                }
                fe
            })
            .collect::<Vec<E::Fr>>()
    }
}
