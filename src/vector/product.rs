use crate::vector::Vector;
use crate::vector::Scalar;
use num_complex::Complex;
use super::{Cartessian, CartessianND};

pub trait Dot<Rhs = Self> : Vector{
    type Output;

    fn dot(&self, rhs : Rhs) -> Self::Output;
}


impl<T : Scalar, const N : usize> Dot<Cartessian<T, N>> for Cartessian<T, N>{
    type Output = T;

    fn dot(&self, rhs : Self) -> T{
        self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T : Scalar, const N : usize> Dot<&'a Cartessian<T, N>> for Cartessian<T, N>{
    type Output = T;

    fn dot(&self, rhs : &'a Cartessian<T, N>) -> T{
        self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T : Scalar, const N : usize> Dot<Cartessian<T, N>> for &'a Cartessian<T, N>{
    type Output = T;

    fn dot(&self, rhs : Cartessian<T, N>) -> T{
        self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T : Scalar, const N : usize> Dot<&'a Cartessian<T, N>> for &'a Cartessian<T, N>{
    type Output = T;

    fn dot(&self, rhs : &'a Cartessian<T, N>) -> T{
        self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
    }
}


impl<T : Scalar> Dot<CartessianND<T>> for CartessianND<T>{
    type Output = T;

    fn dot(&self, rhs : CartessianND<T>) -> T{
        if self.dim() != rhs.dim(){
            panic!("Dot between vectors with different lengths are not available.")
        }
        self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T : Scalar> Dot<&'a CartessianND<T>> for CartessianND<T>{
    type Output = T;

    fn dot(&self, rhs : &'a CartessianND<T>) -> T{
        if self.dim() != rhs.dim(){
            panic!("Dot between vectors with different lengths are not available.")
        }
        self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T : Scalar> Dot<CartessianND<T>> for &'a CartessianND<T>{
    type Output = T;

    fn dot(&self, rhs : CartessianND<T>) -> T{
        if self.dim() != rhs.dim(){
            panic!("Dot between vectors with different lengths are not available.")
        }
        self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T : Scalar> Dot<&'a CartessianND<T>> for &'a CartessianND<T>{
    type Output = T;

    fn dot(&self, rhs : &'a CartessianND<T>) -> T{
        if self.dim() != rhs.dim(){
            panic!("Dot between vectors with different lengths are not available.")
        }
        self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
    }
}



pub trait InnerProduct<RHS=Self>{
    type Output;

    fn inner_product(&self, rhs : RHS) -> Self::Output;
}

macro_rules! impl_inner_product_real {
    ($ty : ident) => {
        impl<const N : usize> InnerProduct<Cartessian<$ty, N>> for Cartessian<$ty, N>{
            type Output = $ty;

            fn inner_product(&self, rhs : Cartessian<$ty, N>) -> Self::Output{
                self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a, const N : usize> InnerProduct<&'a Cartessian<$ty, N>> for Cartessian<$ty, N>{
            type Output = $ty;

            fn inner_product(&self, rhs : &'a Cartessian<$ty, N>) -> Self::Output{
                self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a, const N : usize> InnerProduct<Cartessian<$ty, N>> for &'a Cartessian<$ty, N>{
            type Output = $ty;

            fn inner_product(&self, rhs : Cartessian<$ty, N>) -> Self::Output{
                self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a, const N : usize> InnerProduct<&'a Cartessian<$ty, N>> for &'a Cartessian<$ty, N>{
            type Output = $ty;

            fn inner_product(&self, rhs : &'a Cartessian<$ty, N>) -> Self::Output{
                self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
            }
        }


        impl InnerProduct<CartessianND<$ty>> for CartessianND<$ty>{
            type Output = $ty;

            fn inner_product(&self, rhs : CartessianND<$ty>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a> InnerProduct<&'a CartessianND<$ty>> for CartessianND<$ty>{
            type Output = $ty;

            fn inner_product(&self, rhs : &'a CartessianND<$ty>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a> InnerProduct<CartessianND<$ty>> for &'a CartessianND<$ty>{
            type Output = $ty;

            fn inner_product(&self, rhs : CartessianND<$ty>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a> InnerProduct<&'a CartessianND<$ty>> for &'a CartessianND<$ty>{
            type Output = $ty;

            fn inner_product(&self, rhs : &'a CartessianND<$ty>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
            }
        }

    };
}

impl_inner_product_real!(i8);
impl_inner_product_real!(i16);
impl_inner_product_real!(i32);
impl_inner_product_real!(i64);
impl_inner_product_real!(i128);

impl_inner_product_real!(u8);
impl_inner_product_real!(u16);
impl_inner_product_real!(u32);
impl_inner_product_real!(u64);
impl_inner_product_real!(u128);

impl_inner_product_real!(f32);
impl_inner_product_real!(f64);

macro_rules! impl_inner_product_complex {
    ($ty : ident) => {
        impl<const N : usize> InnerProduct<Cartessian<Complex<$ty>, N>> for Cartessian<Complex<$ty>, N>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : Cartessian<Complex<$ty>, N>) -> Self::Output{
                self.into_iter().zip(rhs.iter()).map(|(x, y)| x.conj() * *y).sum()
            }
        }

        impl<'a, const N : usize> InnerProduct<&'a Cartessian<Complex<$ty>, N>> for Cartessian<Complex<$ty>, N>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : &'a Cartessian<Complex<$ty>, N>) -> Self::Output{
                self.into_iter().zip(rhs).map(|(x, y)| x.conj() * *y).sum()
            }
        }

        impl<'a, const N : usize> InnerProduct<Cartessian<Complex<$ty>, N>> for &'a Cartessian<Complex<$ty>, N>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : Cartessian<Complex<$ty>, N>) -> Self::Output{
                self.into_iter().zip(rhs.iter()).map(|(x, y)| x.conj() * *y).sum()
            }
        }

        impl<'a, const N : usize> InnerProduct<&'a Cartessian<Complex<$ty>, N>> for &'a Cartessian<Complex<$ty>, N>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : &'a Cartessian<Complex<$ty>, N>) -> Self::Output{
                self.into_iter().zip(rhs).map(|(x, y)| x.conj() * *y).sum()
            }
        }


        impl InnerProduct<CartessianND<Complex<$ty>>> for CartessianND<Complex<$ty>>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : CartessianND<Complex<$ty>>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs.iter()).map(|(x, y)| x.conj() * *y).sum()
            }
        }

        impl<'a> InnerProduct<&'a CartessianND<Complex<$ty>>> for CartessianND<Complex<$ty>>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : &'a CartessianND<Complex<$ty>>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs).map(|(x, y)| x.conj() * *y).sum()
            }
        }

        impl<'a> InnerProduct<CartessianND<Complex<$ty>>> for &'a CartessianND<Complex<$ty>>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : CartessianND<Complex<$ty>>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs.iter()).map(|(x, y)| x.conj() * *y).sum()
            }
        }

        impl<'a> InnerProduct<&'a CartessianND<Complex<$ty>>> for &'a CartessianND<Complex<$ty>>{
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs : &'a CartessianND<Complex<$ty>>) -> Self::Output{
                if self.dim() != rhs.dim(){
                    panic!("Inner Product between vectors with different lengths are not available.")
                }
                self.into_iter().zip(rhs).map(|(x, y)| x.conj() * *y).sum()
            }
        }

    };
}

impl_inner_product_complex!(f32);
impl_inner_product_complex!(f64);


pub trait Norm {
    type Output;

    fn norm_l1(&self) -> Self::Output;
    fn norm_l2(&self) -> Self::Output;
}

macro_rules! impl_norm {
    ($ty : ident $(, $tys : ident)*) => {
        $(
        impl<const N : usize> Norm for Cartessian<$tys, N>{
            type Output = $ty;

            #[allow(dead_code)]
            fn norm_l1(&self) -> $ty{
                (self.iter().map(|x| (*x as $ty).abs()).sum::<$ty>())
            }

            #[allow(dead_code)]
            fn norm_l2(&self) -> $ty{
                (self.iter().map(|x| x * x).sum::<$tys>() as $ty).sqrt()
            }
        }

        impl Norm for CartessianND<$tys>{
            type Output = $ty;

            #[allow(dead_code)]
            fn norm_l1(&self) -> $ty{
                (self.iter().map(|x| (*x as $ty).abs()).sum::<$ty>())
            }

            #[allow(dead_code)]
            fn norm_l2(&self) -> $ty{
                (self.iter().map(|x| x * x).sum::<$tys>() as $ty).sqrt()
            }
        }
        )*
    };
}

impl_norm!(f32, f32);
impl_norm!(f64, i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize, f64);

macro_rules! impl_norm_complex {
    ($ty : ident) => {
        impl<const N : usize> Norm for Cartessian<Complex<$ty>, N>{
            type Output = $ty;

            #[allow(dead_code)]
            fn norm_l1(&self) -> $ty{
                self.iter().map(|x| x.norm()).sum::<$ty>()
            }

            #[allow(dead_code)]
            fn norm_l2(&self) -> $ty{
                self.iter().map(|x| x.norm_sqr()).sum::<$ty>().sqrt()
            }
        }

        impl Norm for CartessianND<Complex<$ty>>{
            type Output = $ty;

            #[allow(dead_code)]
            fn norm_l1(&self) -> $ty{
                self.iter().map(|x| x.norm()).sum::<$ty>()
            }

            #[allow(dead_code)]
            fn norm_l2(&self) -> $ty{
                self.iter().map(|x| x.norm_sqr()).sum::<$ty>().sqrt()
            }
        }
    };
}

impl_norm_complex!(f32);
impl_norm_complex!(f64);


#[cfg(test)]
mod test {
    use super::*;
    use num_complex::Complex64;
    use approx::assert_abs_diff_eq;
    use crate::vector::{Cartessian2D, CartessianND};

    #[test]
    fn test_norm(){
        let a = Cartessian2D::from_vec(vec![1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(a.norm_l1(), 1.0);
        assert_abs_diff_eq!(a.norm_l2(), 1.0);

        let b = Cartessian2D::from_vec(vec![1f64, 2f64]).unwrap();
        assert_abs_diff_eq!(b.norm_l1(), 3f64);
        assert_abs_diff_eq!(b.norm_l2(), 5f64.sqrt());

        let c = Cartessian2D::from_vec(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]).unwrap();
        assert_abs_diff_eq!(c.norm_l1(), 5f64.sqrt() + 5f64);
        assert_abs_diff_eq!(c.norm_l2(), 30f64.sqrt());

    }

    #[test]
    fn test_dot_inner_product(){
        let a = Cartessian2D::from_vec(vec![1f64, 0f64]).unwrap();
        let b = Cartessian2D::from_vec(vec![1f64, 2f64]).unwrap();
        assert_eq!(a.dot(&b), 1.0);
        assert_eq!(a.inner_product(&b), 1.0);


        let c = Cartessian2D::from_vec(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]).unwrap();
        let d = Cartessian2D::from_vec(vec![Complex64::new(0.0, 2.0), Complex64::new(3.0, 0.0)]).unwrap();
        assert_eq!(c.dot(&d), Complex64::new(5.0, 14.0));
        assert_eq!(c.inner_product(&d), Complex64::new(13.0, -10.0));
    }


    #[test]
    fn test_norm_nd(){
        let a = CartessianND::new(vec![1f64, 0f64]);
        assert_abs_diff_eq!(a.norm_l2(), 1.0);
        assert_abs_diff_eq!(a.norm_l1(), 1.0);

        let b = CartessianND::new(vec![1f64, 2f64]);
        assert_abs_diff_eq!(b.norm_l1(), 3f64);
        assert_abs_diff_eq!(b.norm_l2(), 5f64.sqrt());

        let c = CartessianND::new(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        assert_abs_diff_eq!(c.norm_l1(), 5f64.sqrt() + 5f64);
        assert_abs_diff_eq!(c.norm_l2(), 30f64.sqrt());
    }

    #[test]
    fn test_dot_inner_product_nd(){
        let a = CartessianND::new(vec![1f64, 0f64]);
        let b = CartessianND::new(vec![1f64, 2f64]);
        assert_eq!(a.dot(&b), 1.0);
        assert_eq!(a.inner_product(&b), 1.0);


        let c = CartessianND::new(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        let d = CartessianND::new(vec![Complex64::new(0.0, 2.0), Complex64::new(3.0, 0.0)]);
        assert_eq!(c.dot(&d), Complex64::new(5.0, 14.0));
        assert_eq!(c.inner_product(&d), Complex64::new(13.0, -10.0));
    }

    #[test]
    #[should_panic]
    fn test_dot_panic(){
        let c = CartessianND::new(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        let d = CartessianND::new(vec![Complex64::new(0.0, 2.0)]);
        c.dot(&d);
    }

    #[test]
    #[should_panic]
    fn test_inner_product_panic(){
        let c = CartessianND::new(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        let d = CartessianND::new(vec![Complex64::new(0.0, 2.0)]);
        c.inner_product(&d);
    }
}
