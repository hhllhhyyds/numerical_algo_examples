use std::ops::{Mul, Range};

use crate::{continuous_func::ContinuousFn, dim1_func::Dim1Fn, fl, float_traits::FromF64};

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
