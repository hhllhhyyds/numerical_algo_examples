use std::ops::{Add, Mul};

use crate::continuous_func::ContinuousFn;

pub trait Dim1Fn<T> {
    fn eval(&self, x: T) -> T;
}

#[derive(Clone, Copy, Debug)]
pub struct Linear<T> {
    pub a: T,
    pub b: T,
}

impl<T> Dim1Fn<T> for Linear<T>
where
    T: Mul<Output = T> + Add<Output = T> + Copy,
{
    fn eval(&self, x: T) -> T {
        self.a * x + self.b
    }
}

impl<T> ContinuousFn for Linear<T> where Linear<T>: Dim1Fn<T> {}
