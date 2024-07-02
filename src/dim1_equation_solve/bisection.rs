use std::ops::{Add, Div, Mul, Sub};

use crate::{
    continuous_func::ContinuousFn,
    dim1_func::Dim1Fn,
    float_traits::{fl, Abs, FloatConst, FromF64, MaxMin},
};

use super::find_root::{FindRootProblem, IterStopCondition, SolveResult, StopReason};

#[allow(clippy::collapsible_else_if)]
pub fn bisection_solve<T, F>(
    problem: &FindRootProblem<T, F>,
    stop_cond: &IterStopCondition<T>,
) -> SolveResult<T>
where
    T: Mul<Output = T>
        + PartialOrd
        + MaxMin
        + Copy
        + Sub<Output = T>
        + Div<Output = T>
        + Mul<Output = T>
        + FloatConst
        + Add<Output = T>
        + Abs
        + FromF64,
    F: ContinuousFn + Dim1Fn<T>,
{
    debug_assert!(problem.is_valid());

    let range = problem.search_range();
    let (mut a, mut b) = (range.start, range.end);

    let mut root = a;
    let mut root_is_a = true;
    let mut root_y = problem.func_eval(root);

    let mut iter_count = 0;

    loop {
        if let Some(limit) = stop_cond.iter_count_limit() {
            if iter_count >= limit {
                return SolveResult {
                    root,
                    iter_count,
                    stop_reason: StopReason::IterCountLimit,
                };
            }
        }

        let c = (a + b) / fl!(2.0);

        if (c - root).abs() <= stop_cond.x_tolorency() * c.abs().max(T::ONE) {
            return SolveResult {
                root: c,
                iter_count,
                stop_reason: StopReason::TolorencyX,
            };
        }

        let f_c = problem.func_eval(c);
        if f_c.abs() <= stop_cond.y_tolorency() {
            return SolveResult {
                root: c,
                iter_count,
                stop_reason: StopReason::TolorencyY,
            };
        }

        if root_y * f_c <= T::ZERO {
            if root_is_a {
                root_is_a = false;
                b = c;
            } else {
                root_is_a = true;
                a = c;
            }
        } else {
            if root_is_a {
                a = c;
            } else {
                b = c;
            }
        }
        root = c;
        root_y = f_c;

        iter_count += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::{AbsDiffEq, RelativeEq};

    use crate::dim1_func::polynomial::Polynomial;

    #[test]
    fn test_meet_tolorency_x() {
        let problem =
            FindRootProblem::new(Polynomial::new(1.0).with_coefficients(&[1.0]), -2.0..3.0);
        let stop_cond = IterStopCondition::new().with_x_tolorency(fl!(1e-6));
        let result = bisection_solve(&problem, &stop_cond);
        assert!(
            result.root().abs_diff_eq(&-1.0, 1e-6),
            "root = {}",
            result.root()
        );
        assert!(result.stop_reason() == StopReason::TolorencyX);
        println!("iter count = {}", result.iter_count());
    }

    #[test]
    fn test_meet_tolorency_y() {
        let f = Polynomial::new(1.0).with_coefficients(&[1.0]);
        let problem = FindRootProblem::new(f.clone(), -2.0..3.0);
        let stop_cond = IterStopCondition::new().with_y_tolorency(fl!(1e-6));
        let result = bisection_solve(&problem, &stop_cond);
        assert!(f.eval(result.root()).abs_diff_eq(&0.0, 1e-6));
        assert!(result.stop_reason() == StopReason::TolorencyY);
        println!("iter count = {}", result.iter_count());
    }

    #[test]
    fn test_meet_iter_count_limit() {
        let f = Polynomial::new(1.0).with_coefficients(&[1.0]);
        let problem = FindRootProblem::new(f, -2.0..3.0);
        let stop_cond = IterStopCondition::new().with_iter_count_limit(10);
        let result = bisection_solve(&problem, &stop_cond);
        assert!(result.stop_reason() == StopReason::IterCountLimit);
        assert!(result.iter_count() == 10);
        assert!(result.root().abs_diff_ne(&-1.0, 1e-6));
    }

    #[test]
    fn test_large_abs_root() {
        let f = Polynomial::new(-1e-4).with_coefficients(&[1.0]);
        let problem = FindRootProblem::new(f, 9000.0..20000.0);
        let stop_cond = IterStopCondition::new().with_x_tolorency(1e-13);
        let result = bisection_solve(&problem, &stop_cond);
        println!("iter count = {}", result.iter_count());
        println!("root = {}", result.root());
        assert!(result.stop_reason() == StopReason::TolorencyX);
        assert!(result.root().relative_eq(&10000.0, 1e-16, 1e-13));
    }
}
