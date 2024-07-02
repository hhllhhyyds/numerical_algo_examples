pub mod polynomial;

pub trait Dim1Fn<T> {
    fn eval(&self, x: T) -> T;
}
