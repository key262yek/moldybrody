use crate::vector::Scalar;
use crate::vector::{
    product::{Dot, InnerProduct, Norm},
    Vector,
};
use crate::vector::{Cartessian, CartessianND};
use approx::AbsDiffEq;
use num_traits::{Float, FloatConst};
use std::convert::From;
use std::fmt::Debug;

use serde::{Deserialize, Serialize};
use std::ops::Rem;
use std::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Polar<T> {
    pub radius: T,
    pub theta: T,
}

impl<T> Vector for Polar<T>
where
    T: Scalar,
{
    type Item = T;

    fn dim(&self) -> usize {
        2
    }
}

impl<T: Scalar> Polar<T> {
    pub fn new(radius: T, theta: T) -> Self
    where
        T: AbsDiffEq<Epsilon = T> + PartialOrd + FloatConst + Float,
    {
        let pi = <T as FloatConst>::PI();
        let (r, t) = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            (T::zero(), T::zero())
        } else if radius < T::zero() {
            (-radius, (theta + pi) % (pi + pi))
        } else {
            (radius, theta % (pi + pi))
        };

        Self {
            radius: r,
            theta: t,
        }
    }

    pub fn get_radius(&self) -> &T {
        &self.radius
    }

    pub fn get_theta(&self) -> &T {
        &self.theta
    }

    pub fn get_mut_radius(&mut self) -> &mut T {
        &mut self.radius
    }

    pub fn get_mut_theta(&mut self) -> &mut T {
        &mut self.theta
    }

    pub fn dim(&self) -> usize {
        2
    }

    pub(crate) fn index(&self, index: usize) -> T
    where
        T: Float,
    {
        match index {
            0 => self.radius * self.theta.cos(),
            1 => self.radius * self.theta.sin(),
            _ => {
                panic!("Out of Bound")
            }
        }
    }
}

impl<T> Default for Polar<T>
where
    T: Default + Scalar,
{
    fn default() -> Self {
        Self {
            radius: T::default(),
            theta: T::default(),
        }
    }
}

impl<T> PartialEq for Polar<T>
where
    T: Scalar + FloatConst + Float,
{
    fn eq(&self, other: &Self) -> bool {
        let dr = self.radius - other.radius;
        let dt = self.theta - other.theta;
        let pi = <T as FloatConst>::PI();

        dr == T::zero() && (dt / (pi + pi)).fract() == T::zero()
    }
}

impl<T> AbsDiffEq for Polar<T>
where
    T: Scalar + AbsDiffEq + FloatConst + Float,
    <T as AbsDiffEq>::Epsilon: Clone,
{
    type Epsilon = <T as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        <T as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let pi = <T as FloatConst>::PI();
        let dt = ((self.theta - other.theta) / pi).fract();
        self.radius.abs_diff_eq(&other.radius, epsilon.clone())
            && dt.abs_diff_eq(&T::zero(), epsilon.clone())
    }
}

impl<T> From<Polar<T>> for Cartessian<T, 2>
where
    T: Float + Scalar,
{
    fn from(polar: Polar<T>) -> Self {
        Self {
            coord: [
                polar.radius * polar.theta.cos(),
                polar.radius * polar.theta.sin(),
            ],
        }
    }
}

impl<T> From<&Polar<T>> for Cartessian<T, 2>
where
    T: Scalar + Float,
{
    fn from(polar: &Polar<T>) -> Self {
        Self {
            coord: [
                polar.radius * polar.theta.cos(),
                polar.radius * polar.theta.sin(),
            ],
        }
    }
}

impl<T> From<&mut Polar<T>> for Cartessian<T, 2>
where
    T: Scalar + Float,
{
    fn from(polar: &mut Polar<T>) -> Self {
        Self {
            coord: [
                polar.radius * polar.theta.cos(),
                polar.radius * polar.theta.sin(),
            ],
        }
    }
}

impl<T> From<Cartessian<T, 2>> for Polar<T>
where
    T: Float + Scalar,
    Cartessian<T, 2>: Norm<Output = T>,
{
    fn from(carte: Cartessian<T, 2>) -> Self {
        let radius = carte.norm_l2();
        let theta = carte[1].atan2(carte[0]);

        Self { radius, theta }
    }
}

impl<T> From<&Cartessian<T, 2>> for Polar<T>
where
    T: Float + Scalar,
    Cartessian<T, 2>: Norm<Output = T>,
{
    fn from(carte: &Cartessian<T, 2>) -> Self {
        let radius = carte.norm_l2();
        let theta = carte[1].atan2(carte[0]);

        Self { radius, theta }
    }
}

impl<T> From<&mut Cartessian<T, 2>> for Polar<T>
where
    T: Float + Scalar,
    Cartessian<T, 2>: Norm<Output = T>,
{
    fn from(carte: &mut Cartessian<T, 2>) -> Self {
        let radius = carte.norm_l2();
        let theta = carte[1].atan2(carte[0]);

        Self { radius, theta }
    }
}

impl<T> From<Polar<T>> for CartessianND<T>
where
    T: Float + Scalar,
{
    fn from(polar: Polar<T>) -> Self {
        Self {
            coord: vec![
                polar.radius * polar.theta.cos(),
                polar.radius * polar.theta.sin(),
            ],
        }
    }
}

impl<T> From<&Polar<T>> for CartessianND<T>
where
    T: Float + Scalar,
{
    fn from(polar: &Polar<T>) -> Self {
        Self {
            coord: vec![
                polar.radius * polar.theta.cos(),
                polar.radius * polar.theta.sin(),
            ],
        }
    }
}

impl<T> From<&mut Polar<T>> for CartessianND<T>
where
    T: Float + Scalar,
{
    fn from(polar: &mut Polar<T>) -> Self {
        Self {
            coord: vec![
                polar.radius * polar.theta.cos(),
                polar.radius * polar.theta.sin(),
            ],
        }
    }
}

impl<T> From<CartessianND<T>> for Polar<T>
where
    T: Float + Scalar,
    CartessianND<T>: Norm<Output = T>,
{
    fn from(carte: CartessianND<T>) -> Self {
        if carte.dim() != 2 {
            panic!("Polar coordinate is avaliable only on 2D domain");
        }
        let radius = carte.norm_l2();
        let theta = carte[1].atan2(carte[0]);

        Self { radius, theta }
    }
}

impl<T> From<&CartessianND<T>> for Polar<T>
where
    T: Float + Scalar,
    CartessianND<T>: Norm<Output = T>,
{
    fn from(carte: &CartessianND<T>) -> Self {
        if carte.dim() != 2 {
            panic!("Polar coordinate is avaliable only on 2D domain");
        }
        let radius = carte.norm_l2();
        let theta = carte[1].atan2(carte[0]);

        Self { radius, theta }
    }
}

impl<T> From<&mut CartessianND<T>> for Polar<T>
where
    T: Float + Scalar,
    CartessianND<T>: Norm<Output = T>,
{
    fn from(carte: &mut CartessianND<T>) -> Self {
        if carte.dim() != 2 {
            panic!("Polar coordinate is avaliable only on 2D domain");
        }
        let radius = carte.norm_l2();
        let theta = carte[1].atan2(carte[0]);

        Self { radius, theta }
    }
}

macro_rules! impl_polar_op {
    ($trt : ident, $operator : tt, $mth : ident) => {
        impl<T> $trt<Polar<T>> for Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Polar<T>;

            fn $mth(mut self, rhs : Polar<T>) -> Polar<T>{
                let x = self.index(0) $operator rhs.index(0);
                let y = self.index(1) $operator rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
                self
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Polar<T>;

            fn $mth(mut self, rhs : &'a Polar<T>) -> Polar<T>{
                let x = self.index(0) $operator rhs.index(0);
                let y = self.index(1) $operator rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
                self
            }
        }

        impl<'a, T> $trt<Polar<T>> for &'a Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Polar<T>;

            fn $mth(self, mut rhs : Polar<T>) -> Polar<T>{
                let x = self.index(0) $operator rhs.index(0);
                let y = self.index(1) $operator rhs.index(1);

                rhs.radius = (x * x + y * y).sqrt();
                rhs.theta = y.atan2(x);
                rhs
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for &'a Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Polar<T>;

            fn $mth(self, rhs : &'a Polar<T>) -> Polar<T>{
                let x = self.index(0) $operator rhs.index(0);
                let y = self.index(1) $operator rhs.index(1);

                Polar{
                    radius : (x * x + y * y).sqrt(),
                    theta : y.atan2(x),
                }
            }
        }


        // ====================================================================

        impl<T> $trt<Cartessian<T, 2>> for Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(self, mut rhs : Cartessian<T, 2>) -> Cartessian<T, 2>{
                for i in 0..2{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a Cartessian<T, 2>> for Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(self, rhs : &'a Cartessian<T, 2>) -> Cartessian<T, 2>{
                let x = self.index(0) $operator *rhs.index(0);
                let y = self.index(1) $operator *rhs.index(1);

                Cartessian{
                    coord : [x, y],
                }
            }
        }

        impl<'a, T> $trt<Cartessian<T, 2>> for &'a Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(self, mut rhs : Cartessian<T, 2>) -> Cartessian<T, 2>{
                for i in 0..2{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a Cartessian<T, 2>> for &'a Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(self, rhs : &'a Cartessian<T, 2>) -> Cartessian<T, 2>{
                let x = self.index(0) $operator *rhs.index(0);
                let y = self.index(1) $operator *rhs.index(1);

                Cartessian{
                    coord : [x, y],
                }
            }
        }

        impl<T> $trt<Polar<T>> for Cartessian<T, 2>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(mut self, rhs : Polar<T>) -> Cartessian<T, 2>{
                for i in 0..2{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for Cartessian<T, 2>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(mut self, rhs : &'a Polar<T>) -> Cartessian<T, 2>{
                for i in 0..2{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<Polar<T>> for &'a Cartessian<T, 2>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(self, rhs : Polar<T>) -> Cartessian<T, 2>{
                let x = *self.index(0) $operator rhs.index(0);
                let y = *self.index(1) $operator rhs.index(1);

                Cartessian{
                    coord : [x, y],
                }
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for &'a Cartessian<T, 2>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = Cartessian<T, 2>;

            fn $mth(self, rhs : &'a Polar<T>) -> Cartessian<T, 2>{
                let x = *self.index(0) $operator rhs.index(0);
                let y = *self.index(1) $operator rhs.index(1);

                Cartessian{
                    coord : [x, y],
                }
            }
        }

        // ====================================================================

        impl<T> $trt<CartessianND<T>> for Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(self, mut rhs : CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                for i in 0..2{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : &'a CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                let x = self.index(0) $operator *rhs.index(0);
                let y = self.index(1) $operator *rhs.index(1);

                CartessianND{
                    coord : vec![x, y],
                }
            }
        }

        impl<'a, T> $trt<CartessianND<T>> for &'a Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(self, mut rhs : CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                for i in 0..2{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for &'a Polar<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : &'a CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                let x = self.index(0) $operator *rhs.index(0);
                let y = self.index(1) $operator *rhs.index(1);

                CartessianND{
                    coord : vec![x, y],
                }
            }
        }

        impl<T> $trt<Polar<T>> for CartessianND<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(mut self, rhs : Polar<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                for i in 0..2{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for CartessianND<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(mut self, rhs : &'a Polar<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                for i in 0..2{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<Polar<T>> for &'a CartessianND<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : Polar<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                let x = *self.index(0) $operator rhs.index(0);
                let y = *self.index(1) $operator rhs.index(1);

                CartessianND{
                    coord : vec![x, y],
                }
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for &'a CartessianND<T>
            where T : $trt<Output = T> + Scalar + Float + Mul<Output =T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : &'a Polar<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                let x = *self.index(0) $operator rhs.index(0);
                let y = *self.index(1) $operator rhs.index(1);

                CartessianND{
                    coord : vec![x, y],
                }
            }
        }

    };
}

impl_polar_op!(Add, +, add);
impl_polar_op!(Sub, -, sub);

impl<T> Mul<T> for Polar<T>
where
    T: Scalar + MulAssign + Clone,
{
    type Output = Polar<T>;

    fn mul(mut self, rhs: T) -> Polar<T> {
        self.radius *= rhs;
        self
    }
}

impl<'a, T> Mul<T> for &'a Polar<T>
where
    T: Scalar + MulAssign + Clone,
{
    type Output = Polar<T>;

    fn mul(self, rhs: T) -> Polar<T> {
        let mut out = self.clone();
        out.radius *= rhs;
        out
    }
}

impl<T> MulAssign<T> for Polar<T>
where
    T: Scalar + MulAssign,
{
    fn mul_assign(&mut self, rhs: T) {
        self.radius *= rhs;
    }
}

impl<T> Div<T> for Polar<T>
where
    T: Scalar + DivAssign + Clone,
{
    type Output = Polar<T>;

    fn div(mut self, rhs: T) -> Polar<T> {
        self.radius /= rhs;
        self
    }
}

impl<'a, T> Div<T> for &'a Polar<T>
where
    T: Scalar + DivAssign + Clone,
{
    type Output = Polar<T>;

    fn div(self, rhs: T) -> Polar<T> {
        let mut out = self.clone();
        out.radius /= rhs;
        out
    }
}

impl<T> DivAssign<T> for Polar<T>
where
    T: Scalar + DivAssign,
{
    fn div_assign(&mut self, rhs: T) {
        self.radius /= rhs;
    }
}

macro_rules! impl_scalar_mul {
    ($ty : ident $(, $tys : ident)*) => {
        impl Mul<Polar<$ty>> for $ty{
            type Output = Polar<$ty>;

            fn mul(self, mut rhs : Polar<$ty>) -> Polar<$ty>{
                rhs.radius *= self;
                rhs
            }
        }

        impl<'a> Mul<&'a Polar<$ty>> for $ty{
            type Output = Polar<$ty>;

            fn mul(self, rhs : &'a Polar<$ty>) -> Polar<$ty>{
                let mut out = rhs.clone();
                out.radius *= self;
                out
            }
        }

        $(
        impl Mul<Polar<$tys>> for $tys{
            type Output = Polar<$tys>;

            fn mul(self, mut rhs : Polar<$tys>) -> Polar<$tys>{
                rhs.radius *= self;
                rhs
            }
        }

        impl<'a> Mul<&'a Polar<$tys>> for $tys{
            type Output = Polar<$tys>;

            fn mul(self, rhs : &'a Polar<$tys>) -> Polar<$tys>{
                let mut out = rhs.clone();
                out.radius *= self;

                out
            }
        }
        )*
    };
}

impl_scalar_mul!(f32, f64);

macro_rules! impl_assign_op{
    ($trt : ident, $operate : tt, $mth : ident, $doc : expr) => {

        impl<T> $trt<Polar<T>> for Polar<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : Polar<T>){
                let x = self.index(0) $operate rhs.index(0);
                let y = self.index(1) $operate rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for Polar<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : &'a Polar<T>){
                let x = self.index(0) $operate rhs.index(0);
                let y = self.index(1) $operate rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
            }
        }

        impl<T> $trt<Cartessian<T, 2>> for Polar<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : Cartessian<T, 2>){
                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
            }
        }

        impl<'a, T> $trt<&'a Cartessian<T, 2>> for Polar<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : &'a Cartessian<T, 2>){
                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
            }
        }

        impl<T> $trt<Polar<T>> for Cartessian<T, 2>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : Polar<T>){
                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for Cartessian<T, 2>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : &'a Polar<T>){
                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
            }
        }

        impl<T> $trt<CartessianND<T>> for Polar<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : CartessianND<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for Polar<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : &'a CartessianND<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);

                self.radius = (x * x + y * y).sqrt();
                self.theta = y.atan2(x);
            }
        }

        impl<T> $trt<Polar<T>> for CartessianND<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : Polar<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
            }
        }

        impl<'a, T> $trt<&'a Polar<T>> for CartessianND<T>
            where T : $trt + Scalar + Float + Mul<Output = T>{
            fn $mth(&mut self, rhs : &'a Polar<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 2 for operate with Polar");
                }

                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
            }
        }
    }
}

impl_assign_op!(AddAssign, +, add_assign, "");
impl_assign_op!(SubAssign, -, sub_assign, "");

impl<T> Neg for Polar<T>
where
    T: Rem<Output = T> + FloatConst + Scalar,
{
    type Output = Self;

    fn neg(mut self) -> Self {
        let pi = T::PI();
        self.theta = (self.theta + pi) % (pi * pi);
        self
    }
}

impl<T> Dot<Polar<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: Polar<T>) -> T {
        self.radius * rhs.radius * (self.theta - rhs.theta).cos()
    }
}

impl<'a, T> Dot<&'a Polar<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: &'a Polar<T>) -> T {
        self.radius * rhs.radius * (self.theta - rhs.theta).cos()
    }
}

impl<T> Dot<Cartessian<T, 2>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: Cartessian<T, 2>) -> T {
        self.index(0) * *rhs.index(0) + self.index(1) * *rhs.index(1)
    }
}

impl<'a, T> Dot<&'a Cartessian<T, 2>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: &'a Cartessian<T, 2>) -> T {
        self.index(0) * *rhs.index(0) + self.index(1) * *rhs.index(1)
    }
}

impl<T> Dot<Polar<T>> for Cartessian<T, 2>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: Polar<T>) -> T {
        *self.index(0) * rhs.index(0) + *self.index(1) * rhs.index(1)
    }
}

impl<'a, T> Dot<&'a Polar<T>> for Cartessian<T, 2>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: &'a Polar<T>) -> T {
        *self.index(0) * rhs.index(0) + *self.index(1) * rhs.index(1)
    }
}

impl<T> Dot<CartessianND<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: CartessianND<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 2 for operate with Polar");
        }

        self.index(0) * *rhs.index(0) + self.index(1) * *rhs.index(1)
    }
}

impl<'a, T> Dot<&'a CartessianND<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: &'a CartessianND<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 2 for operate with Polar");
        }

        self.index(0) * *rhs.index(0) + self.index(1) * *rhs.index(1)
    }
}

impl<T> Dot<Polar<T>> for CartessianND<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: Polar<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 2 for operate with Polar");
        }

        *self.index(0) * rhs.index(0) + *self.index(1) * rhs.index(1)
    }
}

impl<'a, T> Dot<&'a Polar<T>> for CartessianND<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn dot(&self, rhs: &'a Polar<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 2 for operate with Polar");
        }

        *self.index(0) * rhs.index(0) + *self.index(1) * rhs.index(1)
    }
}

impl<T> InnerProduct<Polar<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: Polar<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Polar<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Polar<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<Cartessian<T, 2>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: Cartessian<T, 2>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Cartessian<T, 2>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Cartessian<T, 2>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<Polar<T>> for Cartessian<T, 2>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: Polar<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Polar<T>> for Cartessian<T, 2>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Polar<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<CartessianND<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: CartessianND<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a CartessianND<T>> for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a CartessianND<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<Polar<T>> for CartessianND<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: Polar<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Polar<T>> for CartessianND<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Polar<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> Norm for Polar<T>
where
    T: Scalar + Float,
{
    type Output = T;

    fn norm_l1(&self) -> Self::Output {
        self.radius * (self.theta.cos().abs() + self.theta.sin().abs())
    }

    fn norm_l2(&self) -> Self::Output {
        self.radius
    }

    fn norm_l2_sqr(&self) -> Self::Output {
        self.radius * self.radius
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::{Cartessian, CartessianND};
    use approx::assert_abs_diff_eq;
    use serde_json::{from_str, to_string};
    use std::f64::consts::PI;

    #[test]
    fn test_eq_partialeq_abs_diff_eq() {
        let a = Polar::<f64>::new(1.0f64, PI / 4f64);
        let b = Polar::<f64>::new(1.0f64, PI * 2.25f64);

        assert_eq!(&a, &b);
        a.abs_diff_eq(&b, 1e-6);
    }

    #[test]
    fn test_basic() {
        let mut a = Polar::<f64>::new(1.0f64, PI / 4f64);
        assert_eq!(a.get_radius(), &1.0f64);
        assert_eq!(a.get_theta(), &(PI / 4f64));

        *a.get_mut_radius() = 2.0f64;
        assert_abs_diff_eq!(a.radius, 2.0f64);

        *a.get_mut_theta() = PI / 2f64;
        assert_abs_diff_eq!(a.theta, PI / 2f64);

        assert_eq!(a.dim(), 2);

        let default = Polar::<f64>::default();
        assert_eq!(
            default,
            Polar::<f64> {
                radius: 0f64,
                theta: 0f64
            }
        );

        assert_abs_diff_eq!(a.index(0), 0.0f64);
        assert_abs_diff_eq!(a.index(1), 2.0f64);
    }

    #[test]
    fn test_from() {
        let a = Polar::<f64>::new(1.0f64, PI / 4f64);
        assert_abs_diff_eq!(
            Cartessian::from(&a),
            Cartessian::new([1f64 / 2f64.sqrt(); 2])
        );
        assert_abs_diff_eq!(
            CartessianND::from(&a),
            CartessianND::new(vec![1f64 / 2f64.sqrt(); 2])
        );

        let carte = Cartessian::new([0.0f64, 1.0f64]);
        assert_abs_diff_eq!(
            Polar::<f64>::from(&carte),
            Polar::<f64>::new(1.0f64, PI / 2f64)
        );

        let carte = CartessianND::new(vec![1.0f64, 0.0f64]);
        assert_abs_diff_eq!(
            Polar::<f64>::from(&carte),
            Polar::<f64>::new(1.0f64, 0.0f64)
        );
    }

    #[test]
    #[should_panic]
    fn test_from_panic() {
        let carte = CartessianND::new(vec![1.0f64, 0.0f64, 0.0f64]);
        assert_abs_diff_eq!(
            Polar::<f64>::from(&carte),
            Polar::<f64>::new(1.0f64, 0.0f64)
        );
    }

    #[test]
    fn test_binary_op() {
        let mut a = Polar::<f64>::new(1.0f64, PI / 2f64);
        let b = Polar::<f64>::new(1.0f64, PI);
        let mut carte = Cartessian::new([-1f64, 0f64]);
        let mut carte_nd = CartessianND::new(vec![-1f64, 0f64]);

        (&a + &b).abs_diff_eq(&Polar::<f64>::new(2.0f64.sqrt(), PI * 0.75f64), 1e-6);
        (&a - &b).abs_diff_eq(&Polar::<f64>::new(2.0f64.sqrt(), PI * 0.25f64), 1e-6);
        (&b + &a).abs_diff_eq(&Polar::<f64>::new(2.0f64.sqrt(), PI * 0.75f64), 1e-6);
        (&b - &a).abs_diff_eq(&Polar::<f64>::new(2.0f64.sqrt(), PI * 1.25f64), 1e-6);

        (&a + &carte).abs_diff_eq(&Cartessian::new([-1f64, 1f64]), 1e-6);
        (&a - &carte).abs_diff_eq(&Cartessian::new([1f64, 1f64]), 1e-6);
        (&carte + &a).abs_diff_eq(&Cartessian::new([-1f64, 1f64]), 1e-6);
        (&carte - &a).abs_diff_eq(&Cartessian::new([-1f64, -1f64]), 1e-6);

        (&a + &carte_nd).abs_diff_eq(&CartessianND::new(vec![-1f64, 1f64]), 1e-6);
        (&a - &carte_nd).abs_diff_eq(&CartessianND::new(vec![1f64, 1f64]), 1e-6);
        (&carte_nd + &a).abs_diff_eq(&CartessianND::new(vec![-1f64, 1f64]), 1e-6);
        (&carte_nd - &a).abs_diff_eq(&CartessianND::new(vec![-1f64, -1f64]), 1e-6);

        (&a * 2.0f64).abs_diff_eq(&Polar::<f64>::new(2f64, PI / 2f64), 1e-6);
        (2f64 * &a).abs_diff_eq(&Polar::<f64>::new(2f64, PI / 2f64), 1e-6);
        (&a / 2.0f64).abs_diff_eq(&Polar::<f64>::new(0.5f64, PI / 2f64), 1e-6);

        a += &b;
        a.abs_diff_eq(&Polar::<f64>::new(2.0f64.sqrt(), PI * 0.75f64), 1e-6);
        a -= &b;
        a.abs_diff_eq(&Polar::<f64>::new(1.0f64, PI * 0.5f64), 1e-6);

        a += &carte;
        a.abs_diff_eq(&Polar::<f64>::new(2.0f64.sqrt(), PI * 0.75f64), 1e-6);
        a -= &carte;
        a.abs_diff_eq(&Polar::<f64>::new(1.0f64, PI * 0.5f64), 1e-6);

        a += &carte_nd;
        a.abs_diff_eq(&Polar::<f64>::new(2.0f64.sqrt(), PI * 0.75f64), 1e-6);
        a -= &carte_nd;
        a.abs_diff_eq(&Polar::<f64>::new(1.0f64, PI * 0.5f64), 1e-6);

        carte += a;
        carte.abs_diff_eq(&Cartessian::new([-1f64, 1f64]), 1e-6);
        carte -= a;
        carte.abs_diff_eq(&Cartessian::new([-1f64, 0f64]), 1e-6);

        carte_nd += a;
        carte_nd.abs_diff_eq(&CartessianND::new(vec![-1f64, 1f64]), 1e-6);
        carte_nd -= a;
        carte_nd.abs_diff_eq(&CartessianND::new(vec![-1f64, 0f64]), 1e-6);

        a = -a;
        a.abs_diff_eq(&Polar::<f64>::new(1f64, -PI * 0.5f64), 1e-6);
    }

    #[test]
    #[should_panic]
    fn test_add_panic() {
        let a = Polar::<f64>::new(1.0f64, PI / 2f64);
        let carte = CartessianND::new(vec![-1f64, 0f64, 0f64]);

        let _x = carte + a;
    }

    #[test]
    fn test_dot() {
        let a = Polar::<f64>::new(1.0f64, PI / 2f64);
        let b = Polar::<f64>::new(1.0f64, PI);
        let carte = Cartessian::new([-1f64, 0f64]);
        let carte_nd = CartessianND::new(vec![-1f64, 0f64]);

        assert_abs_diff_eq!(a.dot(&b), 0f64);
        assert_abs_diff_eq!(a.dot(&carte), 0f64);
        assert_abs_diff_eq!(a.dot(&carte_nd), 0f64);
        assert_abs_diff_eq!(carte.dot(&a), 0f64);
        assert_abs_diff_eq!(carte_nd.dot(&a), 0f64);
    }

    #[test]
    #[should_panic]
    fn test_dot_panic() {
        let a = Polar::<f64>::new(1.0f64, PI / 2f64);
        let carte = CartessianND::new(vec![-1f64, 0f64, 0f64]);

        a.dot(&carte);
    }

    #[test]
    fn test_norm() {
        let a = Polar::<f64>::new(1.0f64, PI / 4f64);

        assert_abs_diff_eq!(a.norm_l1(), 2f64.sqrt());
        assert_abs_diff_eq!(a.norm_l2(), 1f64);
    }

    #[test]
    fn test_serde_polar() {
        let a: Polar<f64> = Polar::new(2f64, 0f64);
        let expected = r#"{"radius":2.0,"theta":0.0}"#;
        assert_eq!(expected, to_string(&a).unwrap());

        let expected: Polar<f64> = from_str(&expected).unwrap();
        assert_eq!(a, expected);
    }
}
