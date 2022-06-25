use super::{Cartessian, CartessianND};
use crate::vector::Scalar;
use crate::vector::Vector;
use num_complex::Complex;

pub trait Dot<Rhs = Self>: Vector {
    type Output;

    fn dot(&self, rhs: Rhs) -> Self::Output;
}

impl<T: Scalar, const N: usize> Dot<Cartessian<T, N>> for Cartessian<T, N> {
    type Output = T;

    fn dot(&self, rhs: Self) -> T {
        self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T: Scalar, const N: usize> Dot<&'a Cartessian<T, N>> for Cartessian<T, N> {
    type Output = T;

    fn dot(&self, rhs: &'a Cartessian<T, N>) -> T {
        self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
    }
}

impl<T: Scalar> Dot<CartessianND<T>> for CartessianND<T> {
    type Output = T;

    fn dot(&self, rhs: CartessianND<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dot between vectors with different lengths are not available.")
        }
        self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
    }
}

impl<'a, T: Scalar> Dot<&'a CartessianND<T>> for CartessianND<T> {
    type Output = T;

    fn dot(&self, rhs: &'a CartessianND<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dot between vectors with different lengths are not available.")
        }
        self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
    }
}

pub trait InnerProduct<RHS = Self>: Vector {
    type Output;

    fn inner_product(&self, rhs: RHS) -> Self::Output;
}

macro_rules! impl_inner_product_real {
    ($ty : ident) => {
        impl<const N: usize> InnerProduct<Cartessian<$ty, N>> for Cartessian<$ty, N> {
            type Output = $ty;

            fn inner_product(&self, rhs: Cartessian<$ty, N>) -> Self::Output {
                self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a, const N: usize> InnerProduct<&'a Cartessian<$ty, N>> for Cartessian<$ty, N> {
            type Output = $ty;

            fn inner_product(&self, rhs: &'a Cartessian<$ty, N>) -> Self::Output {
                self.into_iter().zip(rhs).map(|(x, y)| *x * *y).sum()
            }
        }

        impl InnerProduct<CartessianND<$ty>> for CartessianND<$ty> {
            type Output = $ty;

            fn inner_product(&self, rhs: CartessianND<$ty>) -> Self::Output {
                if self.dim() != rhs.dim() {
                    panic!(
                        "Inner Product between vectors with different lengths are not available."
                    )
                }
                self.into_iter().zip(rhs.iter()).map(|(x, y)| *x * *y).sum()
            }
        }

        impl<'a> InnerProduct<&'a CartessianND<$ty>> for CartessianND<$ty> {
            type Output = $ty;

            fn inner_product(&self, rhs: &'a CartessianND<$ty>) -> Self::Output {
                if self.dim() != rhs.dim() {
                    panic!(
                        "Inner Product between vectors with different lengths are not available."
                    )
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
        impl<const N: usize> InnerProduct<Cartessian<Complex<$ty>, N>>
            for Cartessian<Complex<$ty>, N>
        {
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs: Cartessian<Complex<$ty>, N>) -> Self::Output {
                self.into_iter()
                    .zip(rhs.iter())
                    .map(|(x, y)| x.conj() * *y)
                    .sum()
            }
        }

        impl<'a, const N: usize> InnerProduct<&'a Cartessian<Complex<$ty>, N>>
            for Cartessian<Complex<$ty>, N>
        {
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs: &'a Cartessian<Complex<$ty>, N>) -> Self::Output {
                self.into_iter().zip(rhs).map(|(x, y)| x.conj() * *y).sum()
            }
        }

        impl InnerProduct<CartessianND<Complex<$ty>>> for CartessianND<Complex<$ty>> {
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs: CartessianND<Complex<$ty>>) -> Self::Output {
                if self.dim() != rhs.dim() {
                    panic!(
                        "Inner Product between vectors with different lengths are not available."
                    )
                }
                self.into_iter()
                    .zip(rhs.iter())
                    .map(|(x, y)| x.conj() * *y)
                    .sum()
            }
        }

        impl<'a> InnerProduct<&'a CartessianND<Complex<$ty>>> for CartessianND<Complex<$ty>> {
            type Output = Complex<$ty>;

            fn inner_product(&self, rhs: &'a CartessianND<Complex<$ty>>) -> Self::Output {
                if self.dim() != rhs.dim() {
                    panic!(
                        "Inner Product between vectors with different lengths are not available."
                    )
                }
                self.into_iter().zip(rhs).map(|(x, y)| x.conj() * *y).sum()
            }
        }
    };
}

impl_inner_product_complex!(f32);
impl_inner_product_complex!(f64);

pub trait Norm: Vector {
    type Output;

    fn norm_l1(&self) -> Self::Output;
    fn norm_l2(&self) -> Self::Output;
    fn norm_l2_sqr(&self) -> Self::Output;
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

            #[allow(dead_code)]
            fn norm_l2_sqr(&self) -> $ty{
                self.iter().map(|x| x * x).sum::<$tys>() as $ty
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

            #[allow(dead_code)]
            fn norm_l2_sqr(&self) -> $ty{
                self.iter().map(|x| x * x).sum::<$tys>() as $ty
            }
        }
        )*
    };
}

impl_norm!(f32, f32);
impl_norm!(f64, i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize, f64);

macro_rules! impl_norm_complex {
    ($ty : ident) => {
        impl<const N: usize> Norm for Cartessian<Complex<$ty>, N> {
            type Output = $ty;

            #[allow(dead_code)]
            fn norm_l1(&self) -> $ty {
                self.iter().map(|x| x.norm()).sum::<$ty>()
            }

            #[allow(dead_code)]
            fn norm_l2(&self) -> $ty {
                self.iter().map(|x| x.norm_sqr()).sum::<$ty>().sqrt()
            }

            #[allow(dead_code)]
            fn norm_l2_sqr(&self) -> $ty {
                self.iter().map(|x| x.norm_sqr()).sum::<$ty>()
            }
        }

        impl Norm for CartessianND<Complex<$ty>> {
            type Output = $ty;

            #[allow(dead_code)]
            fn norm_l1(&self) -> $ty {
                self.iter().map(|x| x.norm()).sum::<$ty>()
            }

            #[allow(dead_code)]
            fn norm_l2(&self) -> $ty {
                self.iter().map(|x| x.norm_sqr()).sum::<$ty>().sqrt()
            }

            #[allow(dead_code)]
            fn norm_l2_sqr(&self) -> $ty {
                self.iter().map(|x| x.norm_sqr()).sum::<$ty>()
            }
        }
    };
}

impl_norm_complex!(f32);
impl_norm_complex!(f64);

pub trait Cross<Rhs = Self>: Vector {
    type Output: Vector;

    fn cross(&self, rhs: Rhs) -> Self::Output;
}

impl<T> Cross<Cartessian<T, 3>> for Cartessian<T, 3>
where
    T: Scalar + Default,
{
    type Output = Cartessian<T, 3>;

    fn cross(&self, rhs: Cartessian<T, 3>) -> Cartessian<T, 3> {
        let mut res = Cartessian::<T, 3>::default();
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            res[k] += self[i] * rhs[j] - self[j] * rhs[i];
        }
        return res;
    }
}

impl<'a, T> Cross<&'a Cartessian<T, 3>> for Cartessian<T, 3>
where
    T: Scalar + Default,
{
    type Output = Cartessian<T, 3>;

    fn cross(&self, rhs: &'a Cartessian<T, 3>) -> Cartessian<T, 3> {
        let mut res = Cartessian::<T, 3>::default();
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            res[k] += self[i] * rhs[j] - self[j] * rhs[i];
        }
        return res;
    }
}

impl<T> Cross<Cartessian<T, 2>> for Cartessian<T, 2>
where
    T: Scalar + Default,
{
    type Output = Cartessian<T, 1>;

    fn cross(&self, rhs: Cartessian<T, 2>) -> Cartessian<T, 1> {
        return Cartessian::<T, 1>::new([self[0] * rhs[1] - self[1] * rhs[0]]);
    }
}

impl<'a, T> Cross<&'a Cartessian<T, 2>> for Cartessian<T, 2>
where
    T: Scalar + Default,
{
    type Output = Cartessian<T, 1>;

    fn cross(&self, rhs: &'a Cartessian<T, 2>) -> Cartessian<T, 1> {
        return Cartessian::<T, 1>::new([self[0] * rhs[1] - self[1] * rhs[0]]);
    }
}

impl<T> Cross<CartessianND<T>> for CartessianND<T>
where
    T: Scalar + Default,
{
    type Output = CartessianND<T>;

    fn cross(&self, rhs: CartessianND<T>) -> CartessianND<T> {
        if self.dim() != rhs.dim() {
            panic!("Cross product between different length of vector is not available.");
        } else if self.dim() != 3 {
            panic!("Cross product is available only for 3D vector");
        }
        let mut res = CartessianND::<T>::zeros(self.dim());
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            res[k] += self[i] * rhs[j] - self[j] * rhs[i];
        }
        return res;
    }
}

impl<'a, T> Cross<&'a CartessianND<T>> for CartessianND<T>
where
    T: Scalar + Default,
{
    type Output = CartessianND<T>;

    fn cross(&self, rhs: &'a CartessianND<T>) -> CartessianND<T> {
        if self.dim() != rhs.dim() {
            panic!("Cross product between different length of vector is not available.");
        } else if self.dim() != 3 {
            panic!("Cross product is available only for 3D vector");
        }
        let mut res = CartessianND::<T>::zeros(self.dim());
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            res[k] += self[i] * rhs[j] - self[j] * rhs[i];
        }
        return res;
    }
}

pub trait Distance<Rhs = Self> {
    type Output;

    fn distance(self, other: Rhs) -> Self::Output;
}

macro_rules! impl_distance {
    ($ty : ident $(, $tys : ident)*) => {
        $(
        impl<const N : usize> Distance<Cartessian<$tys, N>> for Cartessian<$tys, N>{
            type Output = $ty;

            fn distance(self, other : Cartessian<$tys, N>) -> Self::Output{
                (self.into_iter().zip(other.into_iter()).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }

        impl<'a, const N : usize> Distance<&'a Cartessian<$tys, N>> for Cartessian<$tys, N>{
            type Output = $ty;

            fn distance(self, other : &'a Cartessian<$tys, N>) -> Self::Output{
                (self.into_iter().zip(other).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }

        impl<'a, const N : usize> Distance<Cartessian<$tys, N>> for &'a Cartessian<$tys, N>{
            type Output = $ty;

            fn distance(self, other : Cartessian<$tys, N>) -> Self::Output{
                (self.iter().zip(other.into_iter()).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }

        impl<'a, const N : usize> Distance<&'a Cartessian<$tys, N>> for &'a Cartessian<$tys, N>{
            type Output = $ty;

            fn distance(self, other : &'a Cartessian<$tys, N>) -> Self::Output{
                (self.iter().zip(other).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }

        // ==================================

        impl Distance<CartessianND<$tys>> for CartessianND<$tys>{
            type Output = $ty;

            fn distance(self, other : CartessianND<$tys>) -> Self::Output{
                if self.dim() != other.dim(){
                    panic!("Dimension of vectors are not compatible");
                }

                (self.into_iter().zip(other.into_iter()).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }

        impl<'a> Distance<&'a CartessianND<$tys>> for CartessianND<$tys>{
            type Output = $ty;

            fn distance(self, other : &'a CartessianND<$tys>) -> Self::Output{
                if self.dim() != other.dim(){
                    panic!("Dimension of vectors are not compatible");
                }

                (self.into_iter().zip(other).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }

        impl<'a> Distance<CartessianND<$tys>> for &'a CartessianND<$tys>{
            type Output = $ty;

            fn distance(self, other : CartessianND<$tys>) -> Self::Output{
                if self.dim() != other.dim(){
                    panic!("Dimension of vectors are not compatible");
                }

                (self.iter().zip(other.into_iter()).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }

        impl<'a> Distance<&'a CartessianND<$tys>> for &'a CartessianND<$tys>{
            type Output = $ty;

            fn distance(self, other : &'a CartessianND<$tys>) -> Self::Output{
                if self.dim() != other.dim(){
                    panic!("Dimension of vectors are not compatible");
                }

                (self.iter().zip(other).fold(0 as $tys, |acc, (x, y)| acc + (x - y) * (x - y)) as $ty).sqrt()
            }
        }
        )*
    };
}

impl_distance!(f32, f32);
impl_distance!(f64, i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize, f64);

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::{Cartessian2D, Cartessian3D, CartessianND};
    use approx::assert_abs_diff_eq;
    use num_complex::Complex64;

    #[test]
    fn test_norm() {
        let a = Cartessian2D::new([1f64, 0f64]);
        assert_abs_diff_eq!(a.norm_l1(), 1.0);
        assert_abs_diff_eq!(a.norm_l2(), 1.0);

        let b = Cartessian2D::new([1f64, 2f64]);
        assert_abs_diff_eq!(b.norm_l1(), 3f64);
        assert_abs_diff_eq!(b.norm_l2(), 5f64.sqrt());

        let c = Cartessian2D::new([Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        assert_abs_diff_eq!(c.norm_l1(), 5f64.sqrt() + 5f64);
        assert_abs_diff_eq!(c.norm_l2(), 30f64.sqrt());
    }

    #[test]
    fn test_dot_inner_product() {
        let a = Cartessian2D::new([1f64, 0f64]);
        let b = Cartessian2D::new([1f64, 2f64]);
        assert_eq!(a.dot(&b), 1.0);
        assert_eq!(a.inner_product(&b), 1.0);

        let c = Cartessian2D::new([Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        let d = Cartessian2D::new([Complex64::new(0.0, 2.0), Complex64::new(3.0, 0.0)]);
        assert_eq!(c.dot(&d), Complex64::new(5.0, 14.0));
        assert_eq!(c.inner_product(&d), Complex64::new(13.0, -10.0));
    }

    #[test]
    fn test_norm_nd() {
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
    fn test_dot_inner_product_nd() {
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
    fn test_dot_panic() {
        let c = CartessianND::new(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        let d = CartessianND::new(vec![Complex64::new(0.0, 2.0)]);
        c.dot(&d);
    }

    #[test]
    #[should_panic]
    fn test_inner_product_panic() {
        let c = CartessianND::new(vec![Complex64::new(1.0, 2.0), Complex64::new(3.0, 4.0)]);
        let d = CartessianND::new(vec![Complex64::new(0.0, 2.0)]);
        c.inner_product(&d);
    }

    #[test]
    fn test_cross() {
        let ihat = Cartessian3D::new([1, 0, 0]);
        let jhat = Cartessian3D::new([0, 1, 0]);
        let khat = Cartessian3D::new([0, 0, 1]);

        assert_eq!(ihat.cross(&jhat), khat);
        assert_eq!(jhat.cross(&khat), ihat);
        assert_eq!(khat.cross(&ihat), jhat);
        assert_eq!(ihat.cross(&khat), -jhat.clone());
        assert_eq!(jhat.cross(&ihat), -khat.clone());
        assert_eq!(khat.cross(&jhat), -ihat.clone());
    }

    #[test]
    fn test_distance() {
        let a = Cartessian2D::new([1f64, 0f64]);
        let b = Cartessian2D::new([1f64, 2f64]);
        assert_eq!((&a).distance(&b), 2.0);
        assert_eq!((&a).distance(b.clone()), 2.0);
        assert_eq!(a.clone().distance(&b), 2.0);
        assert_eq!(a.distance(b), 2.0);

        let c = CartessianND::new(vec![1f64, 0f64]);
        let d = CartessianND::new(vec![1f64, 2f64]);
        assert_eq!((&c).distance(&d), 2.0);
        assert_eq!((&c).distance(d.clone()), 2.0);
        assert_eq!(c.clone().distance(&d), 2.0);
        assert_eq!(c.distance(d), 2.0);
    }
}
