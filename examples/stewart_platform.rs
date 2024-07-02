use std::f32::consts::PI;

use macroquad::prelude::*;
use numerical_algo_examples::{
    continuous_func::ContinuousFn,
    dim1_equation_solve::{
        brent::brent_solve,
        find_root::{FindRootProblem, IterStopCondition},
    },
    dim1_func::Dim1Fn,
};

const BIAS: Vec2 = vec2(0.1, 0.1);
const MULTI: f32 = 1000.0;

fn draw_line(p0: Vec2, p1: Vec2, color: Color) {
    let p0 = (p0 + BIAS) * MULTI;
    let p1 = (p1 + BIAS) * MULTI;
    macroquad::prelude::draw_line(p0.x, p0.y, p1.x, p1.y, 5.0, color);
}

fn draw_point(p: Vec2, color: Color) {
    let p = (p + BIAS) * MULTI;
    draw_circle(p.x, p.y, 10., color);
}

#[derive(Clone, Copy)]
struct Triangle {
    l2: f32,
    l3: f32,
    gamma: f32,
    pose_xy: Vec2,
    pose_theta: f32,
}

impl Triangle {
    pub fn p1(&self) -> Vec2 {
        self.pose_xy
    }

    pub fn p2(&self) -> Vec2 {
        self.p1()
            + vec2(
                self.l3 * self.pose_theta.cos(),
                self.l3 * self.pose_theta.sin(),
            )
    }

    pub fn p3(&self) -> Vec2 {
        self.p1()
            + vec2(
                self.l2 * (self.pose_theta + self.gamma).cos(),
                self.l2 * (self.pose_theta + self.gamma).sin(),
            )
    }

    pub fn draw(&self) {
        draw_point(self.p3(), RED);
        draw_point(self.p1(), RED);
        draw_point(self.p2(), RED);
        draw_line(self.p1(), self.p2(), BLACK);
        draw_line(self.p2(), self.p3(), BLACK);
        draw_line(self.p3(), self.p1(), BLACK);
    }
}

#[derive(Clone, Copy)]
struct ThetaFunc {
    x1: f32,
    x2: f32,
    y2: f32,
    gamma: f32,
    l2: f32,
    l3: f32,
    p1: f32,
    p2: f32,
    p3: f32,
}

impl ThetaFunc {
    fn a2(&self, theta: f32) -> f32 {
        self.l3 * theta.cos() - self.x1
    }

    fn b2(&self, theta: f32) -> f32 {
        self.l3 * theta.sin()
    }

    fn a3(&self, theta: f32) -> f32 {
        self.l2 * (theta + self.gamma).cos() - self.x2
    }

    fn b3(&self, theta: f32) -> f32 {
        self.l2 * (theta + self.gamma).sin() - self.y2
    }

    pub fn f(&self, theta: f32) -> (f32, f32, f32) {
        let a2 = self.a2(theta);
        let b2 = self.b2(theta);
        let a3 = self.a3(theta);
        let b3 = self.b3(theta);

        let a2_sq = a2 * a2;
        let b2_sq = b2 * b2;
        let a3_sq = a3 * a3;
        let b3_sq = b3 * b3;
        let p1_sq = self.p1 * self.p1;
        let p2_sq = self.p2 * self.p2;
        let p3_sq = self.p3 * self.p3;

        let j1 = p2_sq - p1_sq - a2_sq - b2_sq;
        let j2 = p3_sq - p1_sq - a3_sq - b3_sq;

        let n1 = b3 * j1 - b2 * j2;
        let n2 = -a3 * j1 + a2 * j2;
        let d = 2.0 * (a2 * b3 - b2 * a3);

        let x = n1 / d;
        let y = n2 / d;
        let f = n1 * n1 + n2 * n2 - p1_sq * d * d;
        (x, y, f)
    }
}

impl Dim1Fn<f32> for ThetaFunc {
    fn eval(&self, x: f32) -> f32 {
        self.f(x).2
    }
}
impl ContinuousFn for ThetaFunc {}

#[macroquad::main("StewartPlatform")]
async fn main() {
    let triangle = Triangle {
        l2: 0.18,
        l3: 0.15,
        gamma: 30_f32.to_radians(),
        pose_xy: vec2(0.2, 0.1),
        pose_theta: 40_f32.to_radians(),
    };

    let fix1 = vec2(0.0, 0.0);
    let fix2 = vec2(0.45, 0.0);
    let fix3 = vec2(0.08, 0.38);

    let len1 = (fix1 - triangle.p1()).length();
    let len2 = (fix2 - triangle.p2()).length();
    let len3 = (fix3 - triangle.p3()).length();

    let mut theta_func = ThetaFunc {
        x1: fix2.x,
        x2: fix3.x,
        y2: fix3.y,
        gamma: triangle.gamma,
        l2: triangle.l2,
        l3: triangle.l3,
        p1: len1,
        p2: len2,
        p3: len3,
    };

    let start = std::time::Instant::now();
    let mut iter_count = 0;
    let mut solve_time = 0.0;
    let mut solve_iter_count = 0;

    loop {
        iter_count += 1;
        clear_background(ORANGE);

        let du = std::time::Instant::now() - start;
        let len1 = len1 * (1.0 + 0.2 * du.as_secs_f32().sin());
        theta_func.p1 = len1;

        let (start_x, start_y) = (-PI, theta_func.eval(-PI));
        assert!(start_y != 0.0);
        let rand_x = loop {
            let x = rand::gen_range(-PI, PI);
            let y = theta_func.eval(x);

            if y * start_y < 0.0 && x != start_x {
                break x;
            }
        };

        let problem = FindRootProblem::new(theta_func, (start_x.min(rand_x))..start_x.max(rand_x));
        let solve_start = std::time::Instant::now();
        let result = brent_solve(
            &problem,
            &IterStopCondition::default()
                .with_x_tolorency(1e-3)
                .with_iter_count_limit(10),
        );
        solve_time += (std::time::Instant::now() - solve_start).as_secs_f64();
        solve_iter_count += result.iter_count();
        if iter_count >= 100 {
            iter_count = 0;
            println!(
                "100 average time used in solve = {:.2} ms, solve iter count = {}",
                solve_time / 100.0 * 1e6,
                solve_iter_count as f32 / 100.0
            );
            solve_time = 0.0;
            solve_iter_count = 0;
        }

        let (x, y, _) = theta_func.f(result.root());
        let mut tri = triangle;
        tri.pose_theta = result.root();
        tri.pose_xy = vec2(x, y);

        tri.draw();

        draw_point(fix1, BLACK);
        draw_point(fix2, BLACK);
        draw_point(fix3, BLACK);

        draw_line(fix1, fix1 + (tri.p1() - fix1).normalize() * len1, BLUE);
        draw_line(fix2, fix2 + (tri.p2() - fix2).normalize() * len2, BLUE);
        draw_line(fix3, fix3 + (tri.p3() - fix3).normalize() * len3, BLUE);

        next_frame().await
    }
}
