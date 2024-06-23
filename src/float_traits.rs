pub trait FromF64 {
    fn from_f64(x: f64) -> Self;
}

#[macro_export(local_inner_macros)]
macro_rules! fl {
    ($x: literal) => {
        FromF64::from_f64($x)
    };
}

impl FromF64 for f32 {
    fn from_f64(x: f64) -> Self {
        x as f32
    }
}

impl FromF64 for f64 {
    fn from_f64(x: f64) -> Self {
        x
    }
}

pub trait FloatConst {
    const EPSILON: Self;
}

impl FloatConst for f32 {
    const EPSILON: Self = f32::EPSILON;
}

impl FloatConst for f64 {
    const EPSILON: Self = f64::EPSILON;
}

pub trait Abs {
    fn abs(self) -> Self;
}

impl Abs for f32 {
    fn abs(self) -> Self {
        f32::abs(self)
    }
}

impl Abs for f64 {
    fn abs(self) -> Self {
        f64::abs(self)
    }
}

pub trait MaxMin: PartialOrd + Sized {
    fn max(self, other: Self) -> Self {
        if self >= other {
            self
        } else {
            other
        }
    }

    fn min(self, other: Self) -> Self {
        if self <= other {
            self
        } else {
            other
        }
    }
}

impl MaxMin for f32 {}
impl MaxMin for f64 {}
