use std::{
    fmt::Display,
    ops::{Add, Mul, Sub},
};

use crate::{
    continuous_func::ContinuousFn,
    float_traits::{Abs, FloatConst},
};

use super::Dim1Fn;

#[derive(Clone, Debug)]
pub struct Polynomial<T> {
    /// coefficient of the highest order
    leading_coe: T,
    coefficients: Option<Vec<T>>,
    base_points: Option<Vec<T>>,
}

impl<T> Polynomial<T> {
    pub fn new(leading_coe: T) -> Self {
        Self {
            leading_coe,
            coefficients: None,
            base_points: None,
        }
    }

    pub fn degree(&self) -> usize {
        if let Some(coes) = &self.coefficients {
            coes.len()
        } else if let Some(bps) = &self.base_points {
            bps.len()
        } else {
            0
        }
    }

    pub fn with_leading_coe(mut self, leading_coe: T) -> Self {
        self.leading_coe = leading_coe;
        self
    }
}

impl<T: FloatConst> Default for Polynomial<T> {
    fn default() -> Self {
        Polynomial::new(T::ONE)
    }
}

impl<T: Clone> Polynomial<T> {
    pub fn with_coefficients(mut self, coefficients: &[T]) -> Self {
        if coefficients.is_empty() {
            return self;
        }
        if let Some(bps) = &self.base_points {
            assert!(bps.len() == coefficients.len());
        }
        self.coefficients = Some(coefficients.to_vec());
        self
    }

    pub fn with_base_points(mut self, base_points: &[T]) -> Self {
        if base_points.is_empty() {
            return self;
        }
        if let Some(coes) = &self.coefficients {
            assert!(base_points.len() == coes.len());
        }
        self.base_points = Some(base_points.to_vec());
        self
    }
}

impl<T> Polynomial<T>
where
    T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Copy,
{
    /// Nested multiplication - An effective algorithm for evaluating the value of polynomial at `x`
    pub fn nest_mul(&self, x: T) -> T {
        let mut y: T = self.leading_coe;

        if let Some(b) = &self.base_points {
            if let Some(c) = &self.coefficients {
                for i in (0..b.len()).rev() {
                    y = y * (x - b[i]) + c[i];
                }
            } else {
                for i in (0..b.len()).rev() {
                    y = y * (x - b[i]);
                }
            }
        } else if let Some(c) = &self.coefficients {
            for i in (0..c.len()).rev() {
                y = y * x + c[i];
            }
        }

        y
    }
}

impl<T> Dim1Fn<T> for Polynomial<T>
where
    T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Copy,
{
    fn eval(&self, x: T) -> T {
        self.nest_mul(x)
    }
}

impl<T> ContinuousFn for Polynomial<T> where Polynomial<T>: Dim1Fn<T> {}

impl<T> Display for Polynomial<T>
where
    T: Display + PartialOrd + Abs + FloatConst,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = format!("{}", self.leading_coe);

        for i in (0..self.degree()).rev() {
            let coe = self.coefficients.as_ref().map(|coes| &coes[i]);
            let bp = self.base_points.as_ref().map(|bps| &bps[i]);

            let mut bp_s = "x".to_string();
            if let Some(bp) = bp {
                if bp.ne(&T::ZERO) {
                    bp_s = "(x ".to_string()
                        + if bp.gt(&T::ZERO) { "-" } else { "+" }
                        + " "
                        + &if bp.gt(&T::ZERO) {
                            format!("{}", bp)
                        } else {
                            format!("{}", bp).replace('-', "")
                        }
                        + ")";
                }
            }
            s = bp_s + " * " + &s;
            if let Some(coe) = coe {
                if coe.ne(&T::ZERO) {
                    s = format!("{coe}")
                        + " "
                        + if coe.gt(&T::ZERO) { "+" } else { "-" }
                        + " "
                        + &s;
                }
                if i != 0 && coe.ne(&T::ZERO) {
                    s = "(".to_string() + &s + ")"
                }
            }
        }
        write!(f, "{s}")?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_display() {
        let mut p = Polynomial::new(8.0).with_coefficients(&[1.0, -2.0, 0.0, 4.0]);
        println!("poly 1 = {p}");
        assert!(p.eval(1.0) == 11.);

        p = p.with_base_points(&[5.0, -6.0, 2.0, 7.0]);
        println!("poly 2 = {p}");

        p = Polynomial::default().with_base_points(&[5.0, -6.0, 0.0, 7.0]);
        println!("poly 3 = {p}");
        assert!(p.eval(1.0) == 168.);

        p = Polynomial::default().with_base_points(&[5.1, -6.2, 0.3, 7.4]);
        println!("poly 4 = {p}");
    }

    #[test]
    fn test_nest_mul() {
        let p = Polynomial::new(2.0).with_coefficients(&[-1.0, 5.0, -3.0, 3.0]);
        assert!(p.nest_mul(0.5) == 1.25);
        assert!(p.nest_mul(-2.0) == -15.0);
        assert!(p.nest_mul(-1.0) == -10.0);
        assert!(p.nest_mul(0.0) == -1.0);
        assert!(p.nest_mul(1.0) == 6.0);
        assert!(p.nest_mul(2.0) == 53.0);

        let p = Polynomial::new(-0.5)
            .with_coefficients(&[1.0, 0.5, 0.5])
            .with_base_points(&[0.0, 2.0, 3.0]);
        assert!(p.nest_mul(1.0) == 0.0)
    }
}
