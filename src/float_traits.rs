pub trait FromF64 {
    fn from_f64(x: f64) -> Self;
}

#[macro_export]
macro_rules! fl {
    ($x: literal) => {
        FromF64::from_f64($x)
    };
}
pub(crate) use fl;

impl FromF64 for f32 {
    #[inline]
    fn from_f64(x: f64) -> Self {
        x as f32
    }
}

impl FromF64 for f64 {
    #[inline]
    fn from_f64(x: f64) -> Self {
        x
    }
}

pub trait FloatConst {
    const EPSILON: Self;
    const ZERO: Self;
    const ONE: Self;
}

impl FloatConst for f32 {
    const EPSILON: Self = f32::EPSILON;
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
}

impl FloatConst for f64 {
    const EPSILON: Self = f64::EPSILON;
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
}

pub trait Abs: std::ops::Neg<Output = Self> + FromF64 + PartialOrd + Sized {
    fn abs(self) -> Self {
        if self > fl!(0.0) {
            self
        } else {
            -self
        }
    }
}

impl Abs for f32 {
    #[inline]
    fn abs(self) -> Self {
        f32::abs(self)
    }
}

impl Abs for f64 {
    #[inline]
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

impl MaxMin for f32 {
    fn max(self, other: Self) -> Self {
        f32::max(self, other)
    }
    fn min(self, other: Self) -> Self {
        f32::min(self, other)
    }
}

impl MaxMin for f64 {
    fn max(self, other: Self) -> Self {
        f64::max(self, other)
    }
    fn min(self, other: Self) -> Self {
        f64::min(self, other)
    }
}
