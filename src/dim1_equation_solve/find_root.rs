use std::ops::{Add, Div, Mul, Range};

use crate::{
    continuous_func::ContinuousFn,
    dim1_func::Dim1Fn,
    fl,
    float_traits::{FloatConst, FromF64},
};

pub struct FindRootProblem<T, F> {
    func: F,
    search_range: Range<T>,
}

impl<T, F> FindRootProblem<T, F> {
    pub fn new_unchecked(func: F, search_range: Range<T>) -> Self {
        Self { func, search_range }
    }
}

impl<T, F> FindRootProblem<T, F>
where
    T: Mul<Output = T> + PartialOrd + FromF64 + Copy,
    F: ContinuousFn + Dim1Fn<T>,
{
    pub fn is_valid(&self) -> bool {
        let f_a = self.func.eval(self.search_range.start);
        let f_b = self.func.eval(self.search_range.end);
        f_a * f_b <= fl!(0.0)
    }

    pub fn new(func: F, search_range: Range<T>) -> Self {
        let p = Self::new_unchecked(func, search_range);
        assert!(p.is_valid());
        p
    }

    pub fn search_range(&self) -> Range<T> {
        self.search_range.clone()
    }

    pub fn func_eval(&self, x: T) -> T {
        self.func.eval(x)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IterStopCondition<T> {
    x_tolorency: T,
    y_tolorency: T,
    iter_count_limit: Option<usize>,
}

impl<T: FloatConst + PartialOrd + FromF64 + Copy> IterStopCondition<T> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_x_tolorency(mut self, x_tol: T) -> Self {
        assert!(x_tol >= fl!(0.0));
        self.x_tolorency = x_tol;
        self
    }

    pub fn with_y_tolorency(mut self, y_tol: T) -> Self {
        assert!(y_tol >= fl!(0.0));
        self.y_tolorency = y_tol;
        self
    }

    pub fn with_iter_count_limit(mut self, limit: usize) -> Self {
        self.iter_count_limit = Some(limit);
        self
    }

    pub fn x_tolorency(&self) -> T {
        self.x_tolorency
    }

    pub fn y_tolorency(&self) -> T {
        self.y_tolorency
    }

    pub fn iter_count_limit(&self) -> Option<usize> {
        self.iter_count_limit
    }
}

impl<T: FloatConst> Default for IterStopCondition<T> {
    fn default() -> Self {
        Self {
            x_tolorency: T::EPSILON,
            y_tolorency: T::EPSILON,
            iter_count_limit: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StopReason {
    TolorencyX,
    TolorencyY,
    IterCountLimit,
}

#[derive(Clone, Debug)]
pub struct SolveResult<T> {
    pub(super) root_range: Range<T>,
    pub(super) iter_count: usize,
    pub(super) stop_reason: StopReason,
}

impl<T> SolveResult<T>
where
    T: Clone,
{
    pub fn root_range(&self) -> Range<T> {
        self.root_range.clone()
    }

    pub fn iter_count(&self) -> usize {
        self.iter_count
    }

    pub fn stop_reason(&self) -> StopReason {
        self.stop_reason
    }
}

impl<T> SolveResult<T>
where
    T: Clone + Add<Output = T> + Div<Output = T> + FromF64 + Copy,
{
    pub fn apporx_root(&self) -> T {
        (self.root_range.start + self.root_range.end) / fl!(2.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::dim1_func;

    use super::*;

    #[test]
    fn test_new() {
        let problem = FindRootProblem::new(dim1_func::Linear { a: 1.0, b: 1.0 }, -1.0..3.0);
        assert!(problem.is_valid());

        let problem =
            FindRootProblem::new_unchecked(dim1_func::Linear { a: 1.0, b: 1.0 }, 1.0..2.0);
        assert!(!problem.is_valid());

        assert!(problem.func_eval(2.0) == 3.0);
    }
}
