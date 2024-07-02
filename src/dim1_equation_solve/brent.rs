use std::{
    fmt::Debug,
    ops::{Add, Div, Mul, Sub},
};

use crate::{
    continuous_func::ContinuousFn,
    dim1_func::Dim1Fn,
    float_traits::{fl, Abs, FloatConst, FromF64, MaxMin},
};

use super::find_root::{FindRootProblem, IterStopCondition, SolveResult, StopReason};

pub fn brent_solve<T, F>(
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
        + Abs
        + Debug,
    F: ContinuousFn + Dim1Fn<T>,
{
    debug_assert!(problem.is_valid());

    let gen_next_point = |next_x| {
        let next_y = problem.func_eval(next_x);
        Point {
            x: next_x,
            y: next_y,
            y_abs: next_y.abs(),
        }
    };

    let range = problem.search_range();
    let (a, b) = (range.start, range.end);
    let c = (a + b) / fl!(2.0);
    let dup = a == c || b == c;

    let mut root_state = RootState {
        left: gen_next_point(a),
        right: gen_next_point(b),
        mid: gen_next_point(c),
        dup,
    };

    let mut iter_count = 0;

    loop {
        if let Some(limit) = stop_cond.iter_count_limit() {
            if iter_count >= limit {
                return SolveResult {
                    root: root_state.mid.x,
                    iter_count,
                    stop_reason: StopReason::IterCountLimit,
                };
            }
        }

        if root_state.right.x - root_state.left.x
            <= stop_cond.x_tolorency() * fl!(2.0) * root_state.mid.x.abs().max(T::ZERO)
        {
            return SolveResult {
                root: root_state.mid.x,
                iter_count,
                stop_reason: StopReason::TolorencyX,
            };
        }

        if root_state.mid.y_abs < stop_cond.y_tolorency() {
            return SolveResult {
                root: root_state.mid.x,
                iter_count,
                stop_reason: StopReason::TolorencyY,
            };
        }

        if !root_state.dup {
            let p = gen_next_point(inverse_quad_interpolation(
                (root_state.left.x, root_state.mid.x, root_state.right.x),
                (root_state.left.y, root_state.mid.y, root_state.right.y),
            ));
            if p.y_abs < stop_cond.y_tolorency() {
                return SolveResult {
                    root: p.x,
                    iter_count,
                    stop_reason: StopReason::TolorencyY,
                };
            }
            if root_state.insert(p) {
                iter_count += 1;
                continue;
            }
        }

        {
            let p = gen_next_point(regula_falsi(
                (root_state.left.x, root_state.right.x),
                (root_state.left.y, root_state.right.y),
            ));
            if p.y_abs < stop_cond.y_tolorency() {
                return SolveResult {
                    root: p.x,
                    iter_count,
                    stop_reason: StopReason::TolorencyY,
                };
            }
            if root_state.insert(p) {
                iter_count += 1;
                continue;
            }
        }

        {
            let next_x = (root_state.left.x + root_state.right.x) / fl!(2.0);
            let p = gen_next_point(next_x);
            if p.y_abs < stop_cond.y_tolorency() {
                return SolveResult {
                    root: p.x,
                    iter_count,
                    stop_reason: StopReason::TolorencyY,
                };
            }

            if root_state.left.y * p.y <= T::ZERO {
                root_state.right = p;
            } else {
                root_state.left = p;
            }
            root_state.mid = gen_next_point((root_state.left.x + root_state.right.x) / fl!(2.0));
            root_state.dup =
                root_state.left.x == root_state.mid.x || root_state.right.x == root_state.mid.x;
            iter_count += 1;
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Point<T> {
    x: T,
    y: T,
    y_abs: T,
}

#[derive(Debug)]
struct RootState<T> {
    left: Point<T>,
    right: Point<T>,
    mid: Point<T>,
    dup: bool, // if mid is same as left or right
}

impl<T> RootState<T>
where
    T: PartialOrd
        + Abs
        + Mul<Output = T>
        + Add<Output = T>
        + Sub<Output = T>
        + Div<Output = T>
        + Copy
        + FloatConst
        + FromF64,
{
    #[allow(clippy::collapsible_else_if)]
    pub fn insert(&mut self, p: Point<T>) -> bool {
        debug_assert!(self.left.y * self.right.y <= T::ZERO);

        let left_mid_same_side = self.left.y * self.mid.y > T::ZERO;
        let left_p_same_side = self.left.y * p.y > T::ZERO;
        let p_mid_not_same_side = left_mid_same_side ^ left_p_same_side;
        let p_right_to_mid = p.x > self.mid.x;
        let p_outside_range = p.x < self.left.x || p.x > self.right.x;

        let old_interval = self.right.x - self.left.x;

        let mut replace = false;

        if p_outside_range {
            if p_right_to_mid {
                // mid as new left, right as new mid, p as new right
                if p_mid_not_same_side
                    && (p.x - self.mid.x) * fl!(2.0) <= old_interval
                    && self.left.y_abs >= p.y_abs
                {
                    replace = true;
                    self.left = self.mid;
                    self.mid = self.right;
                    self.right = p;
                }
            } else {
                // mid as new right, left as new mid, p as new left
                if p_mid_not_same_side
                    && (self.mid.x - p.x) * fl!(2.0) <= old_interval
                    && self.right.y_abs >= p.y_abs
                {
                    replace = true;
                    self.right = self.mid;
                    self.mid = self.left;
                    self.left = p;
                }
            }
        } else {
            if left_mid_same_side {
                if left_p_same_side {
                    if p_right_to_mid {
                        // mid as new left, p as new mid
                        if (self.right.x - self.mid.x) * fl!(2.0) <= old_interval
                            && self.left.y_abs >= p.y_abs
                        {
                            replace = true;
                            self.left = self.mid;
                            self.mid = p;
                        }
                    } else {
                        // p as new left
                        if (self.right.x - p.x) * fl!(2.0) <= old_interval
                            && self.left.y_abs >= p.y_abs
                        {
                            replace = true;
                            self.left = p;
                        }
                    }
                } else {
                    if p_right_to_mid {
                        // p as new right
                        if (p.x - self.left.x) * fl!(2.0) <= old_interval
                            && self.right.y_abs >= p.y_abs
                        {
                            replace = true;
                            self.right = p;
                        }
                    }
                }
            } else {
                if left_p_same_side {
                    if !p_right_to_mid {
                        // p as new left
                        if (self.right.x - p.x) * fl!(2.0) <= old_interval
                            && self.left.y_abs >= p.y_abs
                        {
                            replace = true;
                            self.left = p;
                        }
                    }
                } else {
                    if p_right_to_mid {
                        // p as new right
                        if (p.x - self.left.x) * fl!(2.0) <= old_interval
                            && self.right.y_abs >= p.y_abs
                        {
                            replace = true;
                            self.right = p;
                        }
                    } else {
                        // mid as new right, p as new mid
                        if (self.mid.x - self.left.x) * fl!(2.0) <= old_interval
                            && self.right.y_abs >= p.y_abs
                        {
                            replace = true;
                            self.right = self.mid;
                            self.mid = p;
                        }
                    }
                }
            }
        }

        if replace {
            self.dup = self.left.x == self.mid.x || self.right.x == self.mid.x;
        }

        replace
    }
}

#[allow(non_snake_case)]
#[inline]
fn inverse_quad_interpolation<T>((a, b, c): (T, T, T), (A, B, C): (T, T, T)) -> T
where
    T: Sub<Output = T>
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + PartialEq
        + FloatConst
        + Copy
        + PartialOrd,
{
    debug_assert!(A != B);
    debug_assert!(B != C);
    debug_assert!(A != C);
    let q = A / B;
    let r = C / B;
    let s = C / A;
    c - (r * (r - q) * (c - b) + (T::ONE - r) * s * (c - a))
        / ((q - T::ONE) * (r - T::ONE) * (s - T::ONE))
}

#[allow(non_snake_case)]
#[inline]
fn regula_falsi<T>((a, b): (T, T), (A, B): (T, T)) -> T
where
    T: Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Copy + PartialEq,
{
    debug_assert!(A != B);
    (b * A - a * B) / (A - B)
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::AbsDiffEq;

    use crate::dim1_func::polynomial::Polynomial;

    #[test]
    fn test_meet_tolorency_x() {
        let problem = FindRootProblem::new(
            Polynomial::default().with_coefficients(&[2.0, -3.0]),
            1.5..300.0,
        );
        let stop_cond = IterStopCondition::new().with_y_tolorency(1e-14);
        let result = brent_solve(&problem, &stop_cond);
        println!("iter count = {}", result.iter_count());
        assert!(
            result.root().abs_diff_eq(&2.0, 2.0 * 1e-6)
                || problem.func_eval(result.root()) < f64::EPSILON,
            "root = {}",
            result.root()
        );
    }

    #[test]
    fn test2() {
        let problem = FindRootProblem::new(
            Polynomial::default().with_coefficients(&[-1., 1., 0.]),
            0.0..1.0,
        );
        let result = brent_solve(&problem, &IterStopCondition::new());
        println!("iter count = {}", result.iter_count());
        println!("root = {}", result.root());
        println!("f(root) = {}", problem.func_eval(result.root()));
    }
}
