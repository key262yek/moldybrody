
use crate::vector::basic::Map;
use std::iter::Sum;
use num_traits::{One, Zero};
use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, SubAssign, MulAssign, DivAssign};
use std::fmt::{Debug, Display};
use ndarray_linalg::types::Scalar as NDScalar;
use num_complex::{Complex32, Complex64};
use approx::AbsDiffEq;

use super::{Cartessian, CartessianND};

pub trait Integer : 'static + Copy + Eq + PartialEq + PartialOrd + Ord + Zero + One + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> {}
pub trait Scalar : 'static + Debug + Copy + PartialEq + Sum + Zero + One + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + AddAssign + SubAssign + DivAssign + MulAssign{}
pub trait Float : NDScalar + Display + Debug + Scalar + Send + Sync + AbsDiffEq{}

macro_rules! impl_integer {
    ($name : ident $(, $names : ident)*) => {
        impl Integer for $name {}
        $(
            impl Integer for $names {}
        )*
    };
}

impl_integer!(i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize);

macro_rules! impl_scalar {
    ($name : ident $(, $names : ident)*) => {
        impl Scalar for $name {}
        $(
            impl Scalar for $names {}
        )*
    };
}

impl_scalar!(i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize, f32, f64, Complex32, Complex64);

macro_rules! impl_float {
    ($name : ident $(, $names : ident)*) => {
        impl Float for $name {}
        $(
            impl Float for $names {}
        )*
    };
}

impl_float!(f32, f64, Complex32, Complex64);




fn clone_iopf<A: Copy, B: Copy>(f: impl Fn(A, B) -> A) -> impl FnMut(&mut A, &B) {
    move |x, y| *x = f(*x, *y)
}

fn clone_iopf_rev<A: Copy, B: Copy>(f: impl Fn(A, B) -> B) -> impl FnMut(&mut B, &A) {
    move |x, y| *x = f(*y, *x)
}

macro_rules! impl_binary_op{
    ($trt : ident, $operator : tt, $mth : ident, $iop : tt, $doc : expr) => {
        impl<T, const N : usize> $trt<Cartessian<T, N>> for Cartessian<T, N>
        where T : $trt<Output = T> + Copy + Clone{
            type Output = Cartessian<T, N>;

            fn $mth(mut self, rhs : Self) -> Self{
                self.zip_mut_with(&rhs, clone_iopf(T::$mth));
                self
            }
        }

        impl<'a, T, const N : usize> $trt<&'a Cartessian<T, N>> for Cartessian<T, N>
        where T : $trt<Output = T> + Copy + Clone{
            type Output = Cartessian<T, N>;

            fn $mth(mut self, rhs : &Self) -> Self{
                self.zip_mut_with(rhs, clone_iopf(T::$mth));
                self
            }
        }

        impl<'a, T, const N : usize> $trt<Cartessian<T, N>> for &'a Cartessian<T, N>
        where T : $trt<Output = T> + Copy + Clone{
            type Output = Cartessian<T, N>;

            fn $mth(self, mut rhs : Cartessian<T, N>) -> Cartessian<T, N>{
                rhs.zip_mut_with(self, clone_iopf_rev(T::$mth));
                rhs
            }
        }

        impl<'a, T, const N : usize> $trt<&'a Cartessian<T, N>> for &'a Cartessian<T, N>
        where T : $trt<Output = T> + Copy + Clone + Debug{
            type Output = Cartessian<T, N>;

            fn $mth(self, rhs : &Cartessian<T, N>) -> Cartessian<T, N>{
                let mut out = self.clone();
                out.zip_mut_with(rhs, |x, y| *x = *x $operator *y);
                out
            }
        }

        impl<T, S, const N : usize> $trt<S> for Cartessian<T, N>
        where T : $trt<S, Output = T> + Copy + Clone + Debug,
              S : Scalar{
            type Output = Cartessian<T, N>;

            fn $mth(mut self, c : S) -> Cartessian<T, N>{
                self.map_inplace(move |elt|{
                    *elt = *elt $operator c;
                });
                self
            }
        }

        impl<'a, T, S, const N : usize> $trt<S> for &'a Cartessian<T, N>
        where T : $trt<S, Output = T> + Copy + Clone + Debug,
              S : Scalar{
            type Output = Cartessian<T, N>;

            fn $mth(self, c : S) -> Cartessian<T, N>{
                let mut out = self.clone();
                out.map_inplace( move |elt| {
                    *elt = *elt $operator c;
                });
                out
            }
        }

        impl<'a, T, S, const N : usize> $trt<S> for &'a mut Cartessian<T, N>
        where T : $trt<S, Output = T> + Copy + Clone + Debug,
              S : Scalar{
            type Output = Cartessian<T, N>;

            fn $mth(self, c : S) -> Cartessian<T, N>{
                let mut out = self.clone();
                out.map_inplace( move |elt| {
                    *elt = *elt $operator c;
                });
                out
            }
        }


        impl<T> $trt<CartessianND<T>> for CartessianND<T>
        where T : $trt<Output = T> + Copy + Clone{
            type Output = CartessianND<T>;

            fn $mth(mut self, rhs : Self) -> Self{
                if self.dim() != rhs.dim(){
                    panic!("Binary operation between vectors with different lengths are not available.");
                }

                self.zip_mut_with(&rhs, clone_iopf(T::$mth));
                self
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for CartessianND<T>
        where T : $trt<Output = T> + Copy + Clone{
            type Output = CartessianND<T>;

            fn $mth(mut self, rhs : &Self) -> Self{
                if self.dim() != rhs.dim(){
                    panic!("Binary operation between vectors with different lengths are not available.");
                }

                self.zip_mut_with(rhs, clone_iopf(T::$mth));
                self
            }
        }

        impl<'a, T> $trt<CartessianND<T>> for &'a CartessianND<T>
        where T : $trt<Output = T> + Copy + Clone{
            type Output = CartessianND<T>;

            fn $mth(self, mut rhs : CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim(){
                    panic!("Binary operation between vectors with different lengths are not available.");
                }

                rhs.zip_mut_with(self, clone_iopf_rev(T::$mth));
                rhs
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for &'a CartessianND<T>
        where T : $trt<Output = T> + Copy + Clone + Debug{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : &CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim(){
                    panic!("Binary operation between vectors with different lengths are not available.");
                }

                let mut out = self.clone();
                out.zip_mut_with(rhs, |x, y| *x = *x $operator *y);
                out
            }
        }

        impl<T, S> $trt<S> for CartessianND<T>
        where T : $trt<S, Output = T> + Copy + Clone + Debug,
              S : Scalar{
            type Output = CartessianND<T>;

            fn $mth(mut self, c : S) -> CartessianND<T>{
                self.map_inplace(move |elt|{
                    *elt = *elt $operator c;
                });
                self
            }
        }

        impl<'a, T, S> $trt<S> for &'a CartessianND<T>
        where T : $trt<S, Output = T> + Copy + Clone + Debug,
              S : Scalar{
            type Output = CartessianND<T>;

            fn $mth(self, c : S) -> CartessianND<T>{
                let mut out = self.clone();
                out.map_inplace(move |elt| {
                    *elt = *elt $operator c;
                });
                out
            }
        }

        impl<'a, T, S> $trt<S> for &'a mut CartessianND<T>
        where T : $trt<S, Output = T> + Copy + Clone + Debug,
              S : Scalar{
            type Output = CartessianND<T>;

            fn $mth(self, c : S) -> CartessianND<T>{
                let mut out = self.clone();
                out.map_inplace(move |elt| {
                    *elt = *elt $operator c;
                });
                out
            }
        }
    };
}

macro_rules! if_commutative {
    (Commute { $a:expr } or { $b:expr }) => {
        $a
    };
    (Ordered { $a:expr } or { $b:expr }) => {
        $b
    };
}

macro_rules! impl_scalar_lhs_op {
    ($scalar:ty, $commutative:ident, $operator:tt, $trt:ident, $mth:ident, $doc:expr) => (
        impl<const N : usize> $trt<Cartessian<$scalar, N>> for $scalar{
            type Output = Cartessian<$scalar, N>;
            fn $mth(self, #[allow(unused_mut)] mut rhs: Cartessian<$scalar, N>) -> Cartessian<$scalar, N> {
                if_commutative!($commutative {
                    rhs.$mth(self)
                } or {{
                    rhs.map_inplace(move |elt| {
                        *elt = self $operator *elt;
                    });
                    rhs
                }})
            }
        }

        impl<'a, const N : usize> $trt<&'a Cartessian<$scalar, N>> for $scalar{
            type Output = Cartessian<$scalar, N>;
            fn $mth(self, rhs: &Cartessian<$scalar, N>) -> Cartessian<$scalar, N> {
                if_commutative!($commutative {
                    rhs.$mth(self)
                } or {{
                    let mut out = rhs.clone();
                    out.map_inplace(move |elt|{
                        *elt = self $operator *elt;
                    });
                    out
                }})
            }
        }


        impl $trt<CartessianND<$scalar>> for $scalar{
            type Output = CartessianND<$scalar>;
            fn $mth(self, #[allow(unused_mut)] mut rhs: CartessianND<$scalar>) -> CartessianND<$scalar> {
                if_commutative!($commutative {
                    rhs.$mth(self)
                } or {{
                    rhs.map_inplace(move |elt| {
                        *elt = self $operator *elt;
                    });
                    rhs
                }})
            }
        }

        impl<'a> $trt<&'a CartessianND<$scalar>> for $scalar{
            type Output = CartessianND<$scalar>;
            fn $mth(self, rhs: &CartessianND<$scalar>) -> CartessianND<$scalar> {
                if_commutative!($commutative {
                    rhs.$mth(self)
                } or {{
                    let mut out = rhs.clone();
                    out.map_inplace(move |elt| {
                        *elt = self $operator *elt;
                    });
                    out
                }})
            }
        }
    );
}

impl_binary_op!(Add, +, add, +=, "addition");
impl_binary_op!(Sub, -, sub, -=, "subtraction");
impl_binary_op!(Mul, *, mul, *=, "multiplication");
impl_binary_op!(Div, /, div, /=, "division");

macro_rules! all_scalar_ops {
    ($int_scalar:ty) => (
        impl_scalar_lhs_op!($int_scalar, Commute, +, Add, add, "addition");
        impl_scalar_lhs_op!($int_scalar, Ordered, -, Sub, sub, "subtraction");
        impl_scalar_lhs_op!($int_scalar, Commute, *, Mul, mul, "multiplication");
        impl_scalar_lhs_op!($int_scalar, Ordered, /, Div, div, "division");
    );
}

all_scalar_ops!(i8);
all_scalar_ops!(i16);
all_scalar_ops!(i32);
all_scalar_ops!(i64);
all_scalar_ops!(i128);

all_scalar_ops!(u8);
all_scalar_ops!(u16);
all_scalar_ops!(u32);
all_scalar_ops!(u64);
all_scalar_ops!(u128);

all_scalar_ops!(f32);
all_scalar_ops!(f64);
all_scalar_ops!(Complex32);
all_scalar_ops!(Complex64);



macro_rules! impl_assign_op{
    ($trt : ident, $mth : ident, $doc : expr) => {

        impl<T, const N : usize> $trt<Cartessian<T, N>> for Cartessian<T, N>
        where T : $trt<T> + Copy{

            fn $mth(&mut self, rhs : Cartessian<T, N>){
                self.zip_mut_with(&rhs, |x, y| x.$mth(*y));
            }
        }

        impl<'a, T, const N : usize> $trt<&'a Cartessian<T, N>> for Cartessian<T, N>
        where T : $trt<T> + Copy{

            fn $mth(&mut self, rhs : &'a Cartessian<T, N>){
                self.zip_mut_with(rhs, |x, y| x.$mth(*y));
            }
        }


        impl<'a, T, const N : usize> $trt<T> for Cartessian<T, N>
        where T : $trt<T> + Copy{

            fn $mth(&mut self, rhs : T){
                self.map_inplace(move |elt| {
                    elt.$mth(rhs);
                })
            }
        }

        impl<T> $trt<CartessianND<T>> for CartessianND<T>
        where T : $trt<T> + Copy{

            fn $mth(&mut self, rhs : CartessianND<T>){
                if self.dim() != rhs.dim(){
                    panic!("Binary operation between vectors with different lengths are not available.");
                }

                self.zip_mut_with(&rhs, |x, y| x.$mth(*y));
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for CartessianND<T>
        where T : $trt<T> + Copy{

            fn $mth(&mut self, rhs : &'a CartessianND<T>){
                if self.dim() != rhs.dim(){
                    panic!("Binary operation between vectors with different lengths are not available.");
                }

                self.zip_mut_with(rhs, |x, y| x.$mth(*y));
            }
        }

        impl<'a, T> $trt<T> for CartessianND<T>
        where T : $trt<T> + Copy{

            fn $mth(&mut self, rhs : T){
                self.map_inplace(move |elt| {
                    elt.$mth(rhs);
                })
            }
        }
    };
}

impl_assign_op!(
    AddAssign,
    add_assign,
    "Perform `self += rhs` as elementwise addition (in place).\n"
);
impl_assign_op!(
    SubAssign,
    sub_assign,
    "Perform `self -= rhs` as elementwise subtraction (in place).\n"
);
impl_assign_op!(
    MulAssign,
    mul_assign,
    "Perform `self *= rhs` as elementwise multiplication (in place).\n"
);
impl_assign_op!(
    DivAssign,
    div_assign,
    "Perform `self /= rhs` as elementwise division (in place).\n"
);

impl<T, const N : usize> Neg for Cartessian<T, N>
        where T : Neg<Output = T> + Copy + Clone{
    type Output = Self;

    fn neg(mut self) -> Self{
        self.map_inplace( |elt| {
            *elt = -elt.clone();
        });
        self
    }
}


impl<T> Neg for CartessianND<T>
        where T : Neg<Output = T> + Copy + Clone{
    type Output = Self;

    fn neg(mut self) -> Self{
        self.map_inplace( |elt| {
            *elt = -elt.clone();
        });
        self
    }
}

impl<T, const N : usize> AbsDiffEq<Self> for Cartessian<T, N>
    where T : AbsDiffEq<T>,
          <T as AbsDiffEq>::Epsilon : Clone{
    type Epsilon = T::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        !self.iter().zip(other).any(|(x,y)| !x.abs_diff_eq(y, epsilon.clone()))
    }
}

impl<T> AbsDiffEq<Self> for CartessianND<T>
    where T : AbsDiffEq<T>,
          <T as AbsDiffEq>::Epsilon : Clone{
    type Epsilon = T::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        !self.iter().zip(other).any(|(x,y)| !x.abs_diff_eq(y, epsilon.clone()))
    }
}

impl<T, const N : usize> Cartessian<T, N>
    where T : Scalar{

    pub fn zeros() -> Self{
        Self::new([T::zero(); N])
    }

    pub fn ones() -> Self{
        Self::new([T::one(); N])
    }
}

impl<T> CartessianND<T>
    where T : Scalar{

    pub fn zeros(n : usize) -> Self{
        Self::new(vec![T::zero(); n])
    }

    pub fn ones(n : usize) -> Self{
        Self::new(vec![T::one(); n])
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use crate::vector::Cartessian2D;

    #[test]
    fn test_binary_op(){
        let a = Cartessian2D::from_vec(vec![1, 0]).unwrap();
        let b = Cartessian2D::from_vec(vec![1, 2]).unwrap();

        assert_eq!(a.clone() + b.clone(), Cartessian2D::from_vec(vec![2, 2]).unwrap());
        assert_eq!(a.clone() - b.clone(), Cartessian2D::from_vec(vec![0, -2]).unwrap());
        assert_eq!(a.clone() * b.clone(), Cartessian2D::from_vec(vec![1, 0]).unwrap());

        assert_eq!(a.clone() + &b, Cartessian2D::from_vec(vec![2, 2]).unwrap());
        assert_eq!(a.clone() - &b, Cartessian2D::from_vec(vec![0, -2]).unwrap());
        assert_eq!(a.clone() * &b, Cartessian2D::from_vec(vec![1, 0]).unwrap());

        assert_eq!(&a + b.clone(), Cartessian2D::from_vec(vec![2, 2]).unwrap());
        assert_eq!(&a - b.clone(), Cartessian2D::from_vec(vec![0, -2]).unwrap());
        assert_eq!(&a * b.clone(), Cartessian2D::from_vec(vec![1, 0]).unwrap());

        assert_eq!(&a + &b, Cartessian2D::from_vec(vec![2, 2]).unwrap());
        assert_eq!(&a - &b, Cartessian2D::from_vec(vec![0, -2]).unwrap());
        assert_eq!(&a * &b, Cartessian2D::from_vec(vec![1, 0]).unwrap());

        let a = Cartessian2D::from_vec(vec![1.0, 0.0]).unwrap();
        let b = Cartessian2D::from_vec(vec![1.0, 2.0]).unwrap();

        assert_eq!(a.clone() / b.clone(), Cartessian2D::from_vec(vec![1.0, 0.0]).unwrap());
        assert_eq!(a.clone() / &b, Cartessian2D::from_vec(vec![1.0, 0.0]).unwrap());
        assert_eq!(&a / b.clone(), Cartessian2D::from_vec(vec![1.0, 0.0]).unwrap());
        assert_eq!(&a / &b, Cartessian2D::from_vec(vec![1.0, 0.0]).unwrap());
    }

    #[test]
    fn test_assign_op(){
        let mut a = Cartessian2D::from_vec(vec![1.0, 0.0]).unwrap();
        let b = Cartessian2D::from_vec(vec![1.0, 2.0]).unwrap();

        a += &b;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![2.0, 2.0]).unwrap());

        a += 2f64;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![4.0, 4.0]).unwrap());

        a *= &b;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![4.0, 8.0]).unwrap());

        a *= 2f64;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![8.0, 16.0]).unwrap());

        a -= &b;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![7.0, 14.0]).unwrap());

        a -= 2f64;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![5.0, 12.0]).unwrap());

        a /= &b;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![5.0, 6.0]).unwrap());

        a /= 2f64;
        assert_eq!(&a, &Cartessian2D::from_vec(vec![2.5, 3.0]).unwrap());
    }

    #[test]
    fn test_scalar_op(){
        let a = Cartessian2D::from_vec(vec![1.0, 2.0]).unwrap();
        let c = 3f64;

        assert_abs_diff_eq!(a.clone() + c, Cartessian2D::from_vec(vec![4.0, 5.0]).unwrap());
        assert_abs_diff_eq!(a.clone() - c, Cartessian2D::from_vec(vec![-2.0, -1.0]).unwrap());
        assert_abs_diff_eq!(a.clone() * c, Cartessian2D::from_vec(vec![3.0, 6.0]).unwrap());
        assert_abs_diff_eq!(a.clone() / c, Cartessian2D::from_vec(vec![0.3333333333333333, 0.6666666666666666]).unwrap());

        assert_abs_diff_eq!(&a + c, Cartessian2D::from_vec(vec![4.0, 5.0]).unwrap());
        assert_abs_diff_eq!(&a - c, Cartessian2D::from_vec(vec![-2.0, -1.0]).unwrap());
        assert_abs_diff_eq!(&a * c, Cartessian2D::from_vec(vec![3.0, 6.0]).unwrap());
        assert_abs_diff_eq!(&a / c, Cartessian2D::from_vec(vec![0.3333333333333333, 0.6666666666666666]).unwrap());

        assert_abs_diff_eq!(c + a.clone(), Cartessian2D::from_vec(vec![4.0, 5.0]).unwrap());
        assert_abs_diff_eq!(c - a.clone(), Cartessian2D::from_vec(vec![2.0, 1.0]).unwrap());
        assert_abs_diff_eq!(c * a.clone(), Cartessian2D::from_vec(vec![3.0, 6.0]).unwrap());
        assert_abs_diff_eq!(c / a.clone(), Cartessian2D::from_vec(vec![3.0, 1.5]).unwrap());

        assert_abs_diff_eq!(c + &a, Cartessian2D::from_vec(vec![4.0, 5.0]).unwrap());
        assert_abs_diff_eq!(c - &a, Cartessian2D::from_vec(vec![2.0, 1.0]).unwrap());
        assert_abs_diff_eq!(c * &a, Cartessian2D::from_vec(vec![3.0, 6.0]).unwrap());
        assert_abs_diff_eq!(c / &a, Cartessian2D::from_vec(vec![3.0, 1.5]).unwrap());
    }

    #[test]
    fn test_neg(){
        let a = Cartessian2D::from_vec(vec![1.0, 0.0]).unwrap();
        assert_eq!(-a, Cartessian2D::from_vec(vec![-1.0, 0.0]).unwrap());

        let b = Cartessian2D::from_vec(vec![1.0, 2.0]).unwrap();
        assert_eq!(-b, Cartessian2D::from_vec(vec![-1.0, -2.0]).unwrap());
    }

    #[test]
    fn test_abs_diff_eq(){
        let a = Cartessian2D::from_vec(vec![1f64, 2f64.sqrt()]).unwrap();
        let test = Cartessian2D::from_vec(vec![1f64, 1.4142135623730951]).unwrap();
        assert!(a.abs_diff_eq(&test, 1e-10));
    }

    #[test]
    fn test_binary_op_nd(){
        let a = CartessianND::new(vec![1, 0]);
        let b = CartessianND::new(vec![1, 2]);

        assert_eq!(a.clone() + b.clone(), CartessianND::new(vec![2, 2]));
        assert_eq!(a.clone() - b.clone(), CartessianND::new(vec![0, -2]));
        assert_eq!(a.clone() * b.clone(), CartessianND::new(vec![1, 0]));

        assert_eq!(a.clone() + &b, CartessianND::new(vec![2, 2]));
        assert_eq!(a.clone() - &b, CartessianND::new(vec![0, -2]));
        assert_eq!(a.clone() * &b, CartessianND::new(vec![1, 0]));

        assert_eq!(&a + b.clone(), CartessianND::new(vec![2, 2]));
        assert_eq!(&a - b.clone(), CartessianND::new(vec![0, -2]));
        assert_eq!(&a * b.clone(), CartessianND::new(vec![1, 0]));

        assert_eq!(&a + &b, CartessianND::new(vec![2, 2]));
        assert_eq!(&a - &b, CartessianND::new(vec![0, -2]));
        assert_eq!(&a * &b, CartessianND::new(vec![1, 0]));

        let a = CartessianND::new(vec![1.0, 0.0]);
        let b = CartessianND::new(vec![1.0, 2.0]);

        assert_eq!(a.clone() / b.clone(), CartessianND::new(vec![1.0, 0.0]));
        assert_eq!(a.clone() / &b, CartessianND::new(vec![1.0, 0.0]));
        assert_eq!(&a / b.clone(), CartessianND::new(vec![1.0, 0.0]));
        assert_eq!(&a / &b, CartessianND::new(vec![1.0, 0.0]));
    }

    #[test]
    fn test_assign_op_nd(){
        let mut a = CartessianND::new(vec![1.0, 0.0]);
        let b = CartessianND::new(vec![1.0, 2.0]);

        a += &b;
        assert_eq!(&a, &CartessianND::new(vec![2.0, 2.0]));

        a += 2f64;
        assert_eq!(&a, &CartessianND::new(vec![4.0, 4.0]));

        a *= &b;
        assert_eq!(&a, &CartessianND::new(vec![4.0, 8.0]));

        a *= 2f64;
        assert_eq!(&a, &CartessianND::new(vec![8.0, 16.0]));

        a -= &b;
        assert_eq!(&a, &CartessianND::new(vec![7.0, 14.0]));

        a -= 2f64;
        assert_eq!(&a, &CartessianND::new(vec![5.0, 12.0]));

        a /= &b;
        assert_eq!(&a, &CartessianND::new(vec![5.0, 6.0]));

        a /= 2f64;
        assert_eq!(&a, &CartessianND::new(vec![2.5, 3.0]));
    }

    #[test]
    fn test_scalar_op_nd(){
        let a = CartessianND::new(vec![1.0, 2.0]);
        let c = 3f64;

        assert_abs_diff_eq!(a.clone() + c, CartessianND::new(vec![4.0, 5.0]));
        assert_abs_diff_eq!(a.clone() - c, CartessianND::new(vec![-2.0, -1.0]));
        assert_abs_diff_eq!(a.clone() * c, CartessianND::new(vec![3.0, 6.0]));
        assert_abs_diff_eq!(a.clone() / c, CartessianND::new(vec![0.3333333333333333, 0.6666666666666666]));

        assert_abs_diff_eq!(&a + c, CartessianND::new(vec![4.0, 5.0]));
        assert_abs_diff_eq!(&a - c, CartessianND::new(vec![-2.0, -1.0]));
        assert_abs_diff_eq!(&a * c, CartessianND::new(vec![3.0, 6.0]));
        assert_abs_diff_eq!(&a / c, CartessianND::new(vec![0.3333333333333333, 0.6666666666666666]));

        assert_abs_diff_eq!(c + a.clone(), CartessianND::new(vec![4.0, 5.0]));
        assert_abs_diff_eq!(c - a.clone(), CartessianND::new(vec![2.0, 1.0]));
        assert_abs_diff_eq!(c * a.clone(), CartessianND::new(vec![3.0, 6.0]));
        assert_abs_diff_eq!(c / a.clone(), CartessianND::new(vec![3.0, 1.5]));

        assert_abs_diff_eq!(c + &a, CartessianND::new(vec![4.0, 5.0]));
        assert_abs_diff_eq!(c - &a, CartessianND::new(vec![2.0, 1.0]));
        assert_abs_diff_eq!(c * &a, CartessianND::new(vec![3.0, 6.0]));
        assert_abs_diff_eq!(c / &a, CartessianND::new(vec![3.0, 1.5]));
    }

    #[test]
    fn test_neg_nd(){
        let a = CartessianND::new(vec![1.0, 0.0]);
        assert_eq!(-a, CartessianND::new(vec![-1.0, 0.0]));

        let b = CartessianND::new(vec![1.0, 2.0]);
        assert_eq!(-b, CartessianND::new(vec![-1.0, -2.0]));
    }

    #[test]
    fn test_abs_diff_eq_nd(){
        let a = CartessianND::new(vec![1f64, 2f64.sqrt()]);
        let test = CartessianND::new(vec![1f64, 1.4142135623730951]);
        assert!(a.abs_diff_eq(&test, 1e-10));
    }

    macro_rules! impl_ND_panic {
        ($name : ident, $operator : tt) => {
            #[test]
            #[should_panic]
            fn $name(){
                let a = CartessianND::new(vec![1.0, 0.0]);
                let b = CartessianND::new(vec![1.0, 2.0, 3.0]);

                let _c = a $operator &b;
            }
        };
    }

    impl_ND_panic!(test_panic_add, +);
    impl_ND_panic!(test_panic_sub, -);
    impl_ND_panic!(test_panic_mul, *);
    impl_ND_panic!(test_panic_div, /);

    macro_rules! impl_ND_assign_panic {
        ($name : ident, $operator : tt) => {
            #[test]
            #[should_panic]
            fn $name(){
                let mut a = CartessianND::new(vec![1.0, 0.0]);
                let b = CartessianND::new(vec![1.0, 2.0, 3.0]);

                a $operator &b;
            }
        };
    }

    impl_ND_assign_panic!(test_panic_assign_add, +=);
    impl_ND_assign_panic!(test_panic_assign_sub, -=);
    impl_ND_assign_panic!(test_panic_assign_mul, *=);
    impl_ND_assign_panic!(test_panic_assign_div, /=);
}


