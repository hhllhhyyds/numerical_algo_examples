use std::ops::{Add, Div, Mul, Range, Sub};

use crate::{
    continuous_func::ContinuousFn,
    dim1_func::Dim1Fn,
    fl,
    float_traits::{Abs, FloatConst, FromF64, MaxMin},
};

use super::find_root::FindRootProblem;

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
    root_range: Range<T>,
    iter_count: usize,
    stop_reason: StopReason,
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

pub fn bisection_solve<T, F>(
    problem: &FindRootProblem<T, F>,
    stop_cond: &IterStopCondition<T>,
) -> SolveResult<T>
where
    T: Mul<Output = T>
        + PartialOrd
        + MaxMin
        + FromF64
        + Copy
        + Sub<Output = T>
        + Div<Output = T>
        + Mul<Output = T>
        + FloatConst
        + Add<Output = T>
        + Abs,
    F: ContinuousFn + Dim1Fn<T>,
{
    debug_assert!(problem.is_valid());

    let range = problem.search_range();
    let (mut a, mut b) = (range.start, range.end);

    let mut iter_count = 0;

    let stop_reason = loop {
        if let Some(limit) = stop_cond.iter_count_limit() {
            if iter_count >= limit {
                break StopReason::IterCountLimit;
            }
        }

        if (b - a) <= stop_cond.x_tolorency() * fl!(2.0) * a.max(fl!(1.0)) {
            break StopReason::TolorencyX;
        }

        let c = (a + b) / fl!(2.0);
        let f_c = problem.func_eval(c);
        if f_c.abs() < stop_cond.y_tolorency() {
            break StopReason::TolorencyY;
        }

        let f_a = problem.func_eval(a);
        if f_a * f_c <= fl!(0.0) {
            b = c;
        } else {
            a = c;
        }

        iter_count += 1;
    };

    SolveResult {
        root_range: a..b,
        iter_count,
        stop_reason,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::{AbsDiffEq, RelativeEq};

    use crate::dim1_func;

    #[test]
    fn test_meet_tolorency_x() {
        let problem = FindRootProblem::new(dim1_func::Linear { a: 1.0, b: 1.0 }, -2.0..3.0);
        let stop_cond = IterStopCondition::new().with_x_tolorency(fl!(1e-6));
        let result = bisection_solve(&problem, &stop_cond);
        assert!(
            result.apporx_root().abs_diff_eq(&-1.0, 1e-6),
            "root = {}",
            result.apporx_root()
        );
        assert!(result.stop_reason() == StopReason::TolorencyX);
        println!("iter count = {}", result.iter_count());
    }

    #[test]
    fn test_meet_tolorency_y() {
        let f = dim1_func::Linear { a: 1.0, b: 1.0 };
        let problem = FindRootProblem::new(f, -2.0..3.0);
        let stop_cond = IterStopCondition::new().with_y_tolorency(fl!(1e-6));
        let result = bisection_solve(&problem, &stop_cond);
        assert!(f.eval(result.apporx_root()).abs_diff_eq(&0.0, 1e-6));
        assert!(result.stop_reason() == StopReason::TolorencyY);
        println!("iter count = {}", result.iter_count());
    }

    #[test]
    fn test_meet_iter_count_limit() {
        let f = dim1_func::Linear { a: 1.0, b: 1.0 };
        let problem = FindRootProblem::new(f, -2.0..3.0);
        let stop_cond = IterStopCondition::new().with_iter_count_limit(10);
        let result = bisection_solve(&problem, &stop_cond);
        assert!(result.stop_reason() == StopReason::IterCountLimit);
        assert!(result.iter_count() == 10);
        assert!(result.apporx_root().abs_diff_ne(&-1.0, 1e-6));
    }

    #[test]
    fn test_large_abs_root() {
        let f = dim1_func::Linear { a: -1e-4, b: 1.0 };
        let problem = FindRootProblem::new(f, 9000.0..20000.0);
        let stop_cond = IterStopCondition::new().with_x_tolorency(1e-13);
        let result = bisection_solve(&problem, &stop_cond);
        println!("iter count = {}", result.iter_count());
        println!("root = {}", result.apporx_root());
        assert!(result.stop_reason() == StopReason::TolorencyX);
        assert!(result.apporx_root().relative_eq(&10000.0, 1e-16, 1e-13));
    }
}
