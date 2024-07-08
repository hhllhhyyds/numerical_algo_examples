use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::float_traits::{FloatConst, FromF64};

#[derive(Clone, Copy, Debug)]
pub enum MatrixMajor {
    Row,
    Column,
}

#[derive(Debug, Clone)]
pub struct FullMatrix<T> {
    major: MatrixMajor,
    array: Vec<T>,
    major_dim: usize,
}

impl<T> FullMatrix<T> {
    pub fn new(major: MatrixMajor, major_dim: usize, array: Vec<T>) -> Self {
        assert!(major_dim > 0);
        assert!(!array.is_empty());
        assert!(array.len() % major_dim == 0);
        Self {
            major,
            array,
            major_dim,
        }
    }

    pub fn hilbert(n: usize) -> Self
    where
        T: FloatConst + Div<Output = T> + Add<Output = T> + FromF64,
    {
        let mut v = Vec::default();
        for i in 0..n {
            for j in 0..n {
                v.push(T::ONE / (T::from_f64((i + j) as f64) + T::ONE))
            }
        }
        Self::new(MatrixMajor::Row, n, v)
    }

    fn dim(&self) -> (usize, usize) {
        (self.array.len() / self.major_dim, self.major_dim)
    }

    pub fn row_count(&self) -> usize {
        let dim = self.dim();
        match self.major {
            MatrixMajor::Row => dim.0,
            MatrixMajor::Column => dim.1,
        }
    }

    pub fn column_count(&self) -> usize {
        let dim = self.dim();
        match self.major {
            MatrixMajor::Row => dim.1,
            MatrixMajor::Column => dim.0,
        }
    }

    pub fn is_square_matrix(&self) -> bool {
        self.row_count() == self.column_count()
    }

    pub fn get(&self, non_major: usize, major: usize) -> T
    where
        T: Copy,
    {
        self.array[non_major * self.major_dim + major]
    }

    pub fn set(&mut self, non_major: usize, major: usize, val: T) {
        self.array[non_major * self.major_dim + major] = val;
    }

    pub fn mul_vec(&self, v: &[T]) -> Vec<T>
    where
        T: Copy + Mul<Output = T> + Add<Output = T> + FloatConst,
    {
        match self.major {
            MatrixMajor::Row => {
                assert!(v.len() == self.column_count());
                let mut ret = Vec::default();
                for i in 0..self.row_count() {
                    let mut xi = T::ZERO;
                    for j in 0..self.column_count() {
                        xi = xi + self.get(i, j);
                    }
                    ret.push(xi);
                }

                ret
            }
            MatrixMajor::Column => unimplemented!("not memory efficient for column major matrix"),
        }
    }

    pub fn gaussian_elimination(&self, b: &[T]) -> Vec<T>
    where
        T: Clone
            + Copy
            + Neg<Output = T>
            + Div<Output = T>
            + Mul<Output = T>
            + Add<Output = T>
            + Sub<Output = T>
            + FloatConst,
    {
        match self.major {
            MatrixMajor::Row => {
                assert!(self.is_square_matrix());
                assert!(b.len() == self.column_count());

                let mut b = b.to_owned();
                let mut a = self.clone();
                let n = b.len();

                for col in 0..n {
                    let p = a.get(col, col);
                    for row in (col + 1)..n {
                        let e = -a.get(row, col) / p;
                        a.set(row, col, T::ZERO);
                        for k in (col + 1)..n {
                            let old = a.get(row, k);
                            let sub = a.get(col, k);
                            a.set(row, k, old + e * sub);
                        }
                        b[row] = b[row] + e * b[col];
                    }
                }

                let mut x = vec![T::ZERO; n];
                for i in (0..n).rev() {
                    let mut xi = b[i];
                    for j in ((i + 1)..n).rev() {
                        xi = xi - a.get(i, j) * x[j];
                    }
                    x[i] = xi / a.get(i, i);
                }

                x
            }
            MatrixMajor::Column => unimplemented!("not memory efficient for column major matrix"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gaussian_elimination() {
        let mat = FullMatrix::new(
            MatrixMajor::Row,
            3,
            vec![1., 2., -1., 2., 1., -2., -3., 1., 1.],
        );
        let b = vec![3., 3., -6.];
        let x = mat.gaussian_elimination(&b);
        let anwser = vec![3., 1., 2.];

        for (x0, x1) in x.iter().zip(anwser.iter()) {
            assert!(x0 == x1)
        }
    }

    #[test]
    fn test_gaussian_elimination_solve_hilbert() {
        let hilbert_0 = FullMatrix::hilbert(2);
        println!(
            "solution of hilbert n = 2 is {:?}",
            hilbert_0.gaussian_elimination(&(hilbert_0.mul_vec(&[1.; 2])))
        );

        let hilbert_1 = FullMatrix::hilbert(5);
        println!(
            "solution of hilbert n = 5 is {:?}",
            hilbert_1.gaussian_elimination(&(hilbert_1.mul_vec(&[1.; 5])))
        );

        let hilbert_2 = FullMatrix::hilbert(10);
        println!(
            "solution of hilbert n = 10 is {:?}",
            hilbert_2.gaussian_elimination(&(hilbert_2.mul_vec(&[1.; 10])))
        )
    }
}
