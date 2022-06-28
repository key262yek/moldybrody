use crate::vector::Scalar;
use crate::vector::{
    product::{Dot, InnerProduct, Norm},
    Vector,
};
use crate::vector::{Cartessian, CartessianND};
use approx::AbsDiffEq;
use num_traits::{Float, FloatConst};
use std::convert::From;
use std::convert::TryInto;
use std::fmt::Debug;

use serde::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Rem, Sub, SubAssign};

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Spherical<T> {
    pub radius: T,
    pub theta: T,
    pub phi: T,
}

impl<T> Vector for Spherical<T>
where
    T: Scalar,
{
    type Item = T;

    fn dim(&self) -> usize {
        3
    }
}

fn theta_rep<T>(theta: T) -> T
where
    T: Add<Output = T> + Sub<Output = T> + Rem<Output = T> + PartialOrd + Copy + FloatConst,
{
    let pi = <T as FloatConst>::PI();
    let mut t = theta % (pi + pi);
    if t > pi {
        t = pi + pi - t;
    }
    return t;
}

impl<T> Spherical<T>
where
    T: Scalar,
{
    pub fn new(radius: T, theta: T, phi: T) -> Self
    where
        T: AbsDiffEq<Epsilon = T>
            + Float
            + PartialOrd
            + FloatConst
            + Rem<Output = T>
            + Neg<Output = T>,
    {
        let pi = <T as FloatConst>::PI();
        let (r, t, p) = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            (T::zero(), T::zero(), T::zero())
        } else if radius < T::zero() {
            (-radius, theta_rep(pi - theta), (pi + phi) % (pi + pi))
        } else {
            (radius, theta_rep(theta), phi % (pi + pi))
        };

        Self {
            radius: r,
            theta: t,
            phi: p,
        }
    }

    pub fn get_radius(&self) -> &T {
        &self.radius
    }

    pub fn get_mut_radius(&mut self) -> &mut T {
        &mut self.radius
    }

    pub fn get_theta(&self) -> &T {
        &self.theta
    }

    pub fn get_mut_theta(&mut self) -> &mut T {
        &mut self.theta
    }

    pub fn get_phi(&self) -> &T {
        &self.phi
    }

    pub fn get_mut_phi(&mut self) -> &mut T {
        &mut self.phi
    }

    pub fn dim(&self) -> usize {
        3
    }

    pub(crate) fn index(&self, index: usize) -> T
    where
        T: Float + Scalar,
    {
        match index {
            0 => self.radius * self.theta.sin() * self.phi.cos(),
            1 => self.radius * self.theta.sin() * self.phi.sin(),
            2 => self.radius * self.theta.cos(),
            _ => {
                panic!("Out of Bound")
            }
        }
    }
}

impl<T> Default for Spherical<T>
where
    T: Default + Scalar,
{
    fn default() -> Self {
        Self {
            radius: T::default(),
            theta: T::default(),
            phi: T::default(),
        }
    }
}

impl<T> PartialEq for Spherical<T>
where
    T: PartialEq + Scalar + FloatConst + Rem<Output = T>,
{
    fn eq(&self, other: &Self) -> bool {
        let pi = <T as FloatConst>::PI();
        let dr = self.radius - other.radius;
        let dt = self.theta - other.theta;
        let dt2 = self.theta + other.theta;
        let dp = self.phi - other.phi;

        dr == T::zero()
            && ((dt % (pi + pi) == T::zero()) || (dt2 % (pi + pi) == T::zero()))
            && dp % (pi + pi) == T::zero()
    }
}

impl<T> AbsDiffEq for Spherical<T>
where
    T: AbsDiffEq<Epsilon = T> + Scalar + FloatConst + Rem<Output = T>,
{
    type Epsilon = <T as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        <T as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let pi = <T as FloatConst>::PI();
        let dt = self.theta - other.theta;
        let dt2 = self.theta + other.theta;
        let dp = self.phi - other.phi;

        self.radius.abs_diff_eq(&other.radius, epsilon.clone())
            && ((dt % (pi + pi)).abs_diff_eq(&T::zero(), epsilon.clone())
                || (dt2 % (pi + pi)).abs_diff_eq(&T::zero(), epsilon.clone()))
            && (dp % (pi + pi)).abs_diff_eq(&T::zero(), epsilon.clone())
    }
}

impl<T> From<Spherical<T>> for Cartessian<T, 3>
where
    T: Scalar + Float,
{
    fn from(spheric: Spherical<T>) -> Self {
        Self {
            coord: [spheric.index(0), spheric.index(1), spheric.index(2)],
        }
    }
}

impl<T> From<&Spherical<T>> for Cartessian<T, 3>
where
    T: Scalar + Float,
{
    fn from(spheric: &Spherical<T>) -> Self {
        Self {
            coord: [spheric.index(0), spheric.index(1), spheric.index(2)],
        }
    }
}

impl<T> From<&mut Spherical<T>> for Cartessian<T, 3>
where
    T: Scalar + Float,
{
    fn from(spheric: &mut Spherical<T>) -> Self {
        Self {
            coord: [spheric.index(0), spheric.index(1), spheric.index(2)],
        }
    }
}

impl<T> From<Cartessian<T, 3>> for Spherical<T>
where
    T: AbsDiffEq<Epsilon = T> + Scalar + Float,
    Cartessian<T, 3>: Norm<Output = T>,
{
    fn from(carte: Cartessian<T, 3>) -> Self {
        let radius = carte.norm_l2();
        let theta = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            T::zero()
        } else {
            (carte[2] / radius).acos()
        };
        let phi = carte[1].atan2(carte[0]);

        Self { radius, theta, phi }
    }
}

impl<T> From<&Cartessian<T, 3>> for Spherical<T>
where
    T: AbsDiffEq<Epsilon = T> + Scalar + Float,
    Cartessian<T, 3>: Norm<Output = T>,
{
    fn from(carte: &Cartessian<T, 3>) -> Self {
        let radius = carte.norm_l2();
        let theta = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            T::zero()
        } else {
            (carte[2] / radius).acos()
        };
        let phi = carte[1].atan2(carte[0]);

        Self { radius, theta, phi }
    }
}

impl<T> From<&mut Cartessian<T, 3>> for Spherical<T>
where
    T: AbsDiffEq<Epsilon = T> + Scalar + Float,
    Cartessian<T, 3>: Norm<Output = T>,
{
    fn from(carte: &mut Cartessian<T, 3>) -> Self {
        let radius = carte.norm_l2();
        let theta = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            T::zero()
        } else {
            (carte[2] / radius).acos()
        };
        let phi = carte[1].atan2(carte[0]);

        Self { radius, theta, phi }
    }
}

impl<T> From<Spherical<T>> for CartessianND<T>
where
    T: Scalar + Float,
{
    fn from(spheric: Spherical<T>) -> Self {
        Self {
            coord: vec![spheric.index(0), spheric.index(1), spheric.index(2)],
        }
    }
}

impl<T> From<&Spherical<T>> for CartessianND<T>
where
    T: Scalar + Float,
{
    fn from(spheric: &Spherical<T>) -> Self {
        Self {
            coord: vec![spheric.index(0), spheric.index(1), spheric.index(2)],
        }
    }
}

impl<T> From<&mut Spherical<T>> for CartessianND<T>
where
    T: Scalar + Float,
{
    fn from(spheric: &mut Spherical<T>) -> Self {
        Self {
            coord: vec![spheric.index(0), spheric.index(1), spheric.index(2)],
        }
    }
}

impl<T> From<CartessianND<T>> for Spherical<T>
where
    T: AbsDiffEq<Epsilon = T> + Scalar + Float,
    CartessianND<T>: Norm<Output = T>,
{
    fn from(carte: CartessianND<T>) -> Self {
        if carte.dim() != 3 {
            panic!("Spherical coordinate is avaliable only on 3D domain");
        }
        let radius = carte.norm_l2();
        let theta = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            T::zero()
        } else {
            (carte[2] / radius).acos()
        };
        let phi = carte[1].atan2(carte[0]);

        Self { radius, theta, phi }
    }
}

impl<T> From<&CartessianND<T>> for Spherical<T>
where
    T: AbsDiffEq<Epsilon = T> + Scalar + Float,
    CartessianND<T>: Norm<Output = T>,
{
    fn from(carte: &CartessianND<T>) -> Self {
        if carte.dim() != 3 {
            panic!("Spherical coordinate is avaliable only on 3D domain");
        }
        let radius = carte.norm_l2();
        let theta = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            T::zero()
        } else {
            (carte[2] / radius).acos()
        };
        let phi = carte[1].atan2(carte[0]);

        Self { radius, theta, phi }
    }
}

impl<T> From<&mut CartessianND<T>> for Spherical<T>
where
    T: AbsDiffEq<Epsilon = T> + Scalar + Float,
    CartessianND<T>: Norm<Output = T>,
{
    fn from(carte: &mut CartessianND<T>) -> Self {
        if carte.dim() != 3 {
            panic!("Spherical coordinate is avaliable only on 3D domain");
        }
        let radius = carte.norm_l2();
        let theta = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {
            T::zero()
        } else {
            (carte[2] / radius).acos()
        };
        let phi = carte[1].atan2(carte[0]);

        Self { radius, theta, phi }
    }
}

macro_rules! impl_spherical_op {
    ($trt : ident, $operator : tt, $mth : ident) => {
        impl<T> $trt<Spherical<T>> for Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Spherical<T>;

            fn $mth(mut self, rhs : Spherical<T>) -> Spherical<T>{
                let coord = (0..3).map(|i| self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>();

                self.radius = (coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(coord[2] / self.radius).acos()};
                self.phi = coord[1].atan2(coord[0]);
                self
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Spherical<T>;

            fn $mth(mut self, rhs : &'a Spherical<T>) -> Spherical<T>{
                let coord = (0..3).map(|i| self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>();

                self.radius = (coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(coord[2] / self.radius).acos()};
                self.phi = coord[1].atan2(coord[0]);
                self

            }
        }

        impl<'a, T> $trt<Spherical<T>> for &'a Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Spherical<T>;

            fn $mth(self, mut rhs : Spherical<T>) -> Spherical<T>{
                let coord = (0..3).map(|i| self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>();

                rhs.radius = (coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]).sqrt();
                rhs.theta = if rhs.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(coord[2] / rhs.radius).acos()};
                rhs.phi = coord[1].atan2(coord[0]);
                rhs
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for &'a Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Spherical<T>;

            fn $mth(self, rhs : &'a Spherical<T>) -> Spherical<T>{
                let coord = (0..3).map(|i| self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>();

                let radius = (coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]).sqrt();
                let theta = if radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(coord[2] / radius).acos()};
                let phi = coord[1].atan2(coord[0]);

                Spherical{
                    radius,
                    theta,
                    phi,
                }
            }
        }


        // ====================================================================

        impl<T> $trt<Cartessian<T, 3>> for Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(self, mut rhs : Cartessian<T, 3>) -> Cartessian<T, 3>{
                for i in 0..3{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a Cartessian<T, 3>> for Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(self, rhs : &'a Cartessian<T, 3>) -> Cartessian<T, 3>{
                Cartessian{
                    coord : (0..3).map(|i| self.index(i) $operator *rhs.index(i))
                                  .collect::<Vec<T>>().try_into().unwrap(),
                }
            }
        }

        impl<'a, T> $trt<Cartessian<T, 3>> for &'a Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(self, mut rhs : Cartessian<T, 3>) -> Cartessian<T, 3>{
                for i in 0..3{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a Cartessian<T, 3>> for &'a Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(self, rhs : &'a Cartessian<T, 3>) -> Cartessian<T, 3>{
                Cartessian{
                    coord : (0..3).map(|i| self.index(i) $operator *rhs.index(i))
                                  .collect::<Vec<T>>().try_into().unwrap(),
                }
            }
        }

        impl<T> $trt<Spherical<T>> for Cartessian<T, 3>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(mut self, rhs : Spherical<T>) -> Cartessian<T, 3>{
                for i in 0..3{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for Cartessian<T, 3>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(mut self, rhs : &'a Spherical<T>) -> Cartessian<T, 3>{
                for i in 0..3{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<Spherical<T>> for &'a Cartessian<T, 3>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(self, rhs : Spherical<T>) -> Cartessian<T, 3>{
                Cartessian{
                    coord : (0..3).map(|i| *self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>().try_into().unwrap(),
                }
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for &'a Cartessian<T, 3>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = Cartessian<T, 3>;

            fn $mth(self, rhs : &'a Spherical<T>) -> Cartessian<T, 3>{
                Cartessian{
                    coord : (0..3).map(|i| *self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>().try_into().unwrap(),
                }
            }
        }

        // ====================================================================

        impl<T> $trt<CartessianND<T>> for Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(self, mut rhs : CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                for i in 0..3{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : &'a CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                CartessianND{
                    coord : (0..3).map(|i| self.index(i) $operator *rhs.index(i))
                                  .collect::<Vec<T>>(),
                }
            }
        }

        impl<'a, T> $trt<CartessianND<T>> for &'a Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(self, mut rhs : CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                for i in 0..3{
                    rhs[i] = self.index(i) $operator *rhs.index(i);
                }
                rhs
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for &'a Spherical<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : &'a CartessianND<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                CartessianND{
                    coord : (0..3).map(|i| self.index(i) $operator *rhs.index(i))
                                  .collect::<Vec<T>>(),
                }
            }
        }

        impl<T> $trt<Spherical<T>> for CartessianND<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(mut self, rhs : Spherical<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                for i in 0..3{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for CartessianND<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(mut self, rhs : &'a Spherical<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                for i in 0..3{
                    self[i] = *self.index(i) $operator rhs.index(i);
                }
                self
            }
        }

        impl<'a, T> $trt<Spherical<T>> for &'a CartessianND<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : Spherical<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                CartessianND{
                    coord : (0..3).map(|i| *self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>(),
                }
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for &'a CartessianND<T>
            where T : Scalar + $trt<Output = T> + Float + AbsDiffEq<Epsilon = T>{
            type Output = CartessianND<T>;

            fn $mth(self, rhs : &'a Spherical<T>) -> CartessianND<T>{
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                CartessianND{
                    coord : (0..3).map(|i| *self.index(i) $operator rhs.index(i))
                                  .collect::<Vec<T>>(),
                }
            }
        }

    };
}

impl_spherical_op!(Add, +, add);
impl_spherical_op!(Sub, -, sub);

impl<T> Mul<T> for Spherical<T>
where
    T: Scalar + MulAssign + Clone,
{
    type Output = Spherical<T>;

    fn mul(mut self, rhs: T) -> Spherical<T> {
        self.radius *= rhs;
        self
    }
}

impl<'a, T> Mul<T> for &'a Spherical<T>
where
    T: Scalar + MulAssign + Clone,
{
    type Output = Spherical<T>;

    fn mul(self, rhs: T) -> Spherical<T> {
        let mut out = self.clone();
        out.radius *= rhs;
        out
    }
}

impl<T> MulAssign<T> for Spherical<T>
where
    T: Scalar + MulAssign,
{
    fn mul_assign(&mut self, rhs: T) {
        self.radius *= rhs;
    }
}

impl<T> Div<T> for Spherical<T>
where
    T: Scalar + DivAssign + Clone,
{
    type Output = Spherical<T>;

    fn div(mut self, rhs: T) -> Spherical<T> {
        self.radius /= rhs;
        self
    }
}

impl<'a, T> Div<T> for &'a Spherical<T>
where
    T: Scalar + DivAssign + Clone,
{
    type Output = Spherical<T>;

    fn div(self, rhs: T) -> Spherical<T> {
        let mut out = self.clone();
        out.radius /= rhs;
        out
    }
}

impl<T> DivAssign<T> for Spherical<T>
where
    T: Scalar + DivAssign,
{
    fn div_assign(&mut self, rhs: T) {
        self.radius /= rhs;
    }
}

macro_rules! impl_scalar_op_spherical {
    ($ty : ident) => {
        impl Mul<Spherical<$ty>> for $ty {
            type Output = Spherical<$ty>;

            fn mul(self, mut rhs: Spherical<$ty>) -> Spherical<$ty> {
                rhs.radius *= self;
                rhs
            }
        }

        impl<'a> Mul<&'a Spherical<$ty>> for $ty {
            type Output = Spherical<$ty>;

            fn mul(self, rhs: &'a Spherical<$ty>) -> Spherical<$ty> {
                let mut out = rhs.clone();
                out.radius *= self;
                out
            }
        }
    };
}

impl_scalar_op_spherical!(f32);
impl_scalar_op_spherical!(f64);

macro_rules! impl_assign_op{
    ($trt : ident, $operate : tt, $mth : ident, $doc : expr) => {

        impl<T> $trt<Spherical<T>> for Spherical<T>
            where T : Float + Scalar + AbsDiffEq<Epsilon = T>{
            fn $mth(&mut self, rhs : Spherical<T>){
                let x = self.index(0) $operate rhs.index(0);
                let y = self.index(1) $operate rhs.index(1);
                let z = self.index(2) $operate rhs.index(2);

                self.radius = (x * x + y * y + z * z).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(z / self.radius).acos()};
                self.phi = y.atan2(x);
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for Spherical<T>
            where T : Float + Scalar + AbsDiffEq<Epsilon = T>{
            fn $mth(&mut self, rhs : &'a Spherical<T>){
                let x = self.index(0) $operate rhs.index(0);
                let y = self.index(1) $operate rhs.index(1);
                let z = self.index(2) $operate rhs.index(2);

                self.radius = (x * x + y * y + z * z).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(z / self.radius).acos()};
                self.phi = y.atan2(x);
            }
        }

        impl<T> $trt<Cartessian<T, 3>> for Spherical<T>
            where T : Float + Scalar + AbsDiffEq<Epsilon = T>{
            fn $mth(&mut self, rhs : Cartessian<T, 3>){
                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);
                let z = self.index(2) $operate *rhs.index(2);

                self.radius = (x * x + y * y + z * z).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(z / self.radius).acos()};
                self.phi = y.atan2(x);
            }
        }

        impl<'a, T> $trt<&'a Cartessian<T, 3>> for Spherical<T>
            where T : Float + Scalar + AbsDiffEq<Epsilon = T>{
            fn $mth(&mut self, rhs : &'a Cartessian<T, 3>){
                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);
                let z = self.index(2) $operate *rhs.index(2);

                self.radius = (x * x + y * y + z * z).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(z / self.radius).acos()};
                self.phi = y.atan2(x);
            }
        }

        impl<T> $trt<Spherical<T>> for Cartessian<T, 3>
            where T : Float + Scalar{
            fn $mth(&mut self, rhs : Spherical<T>){
                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
                self[2].$mth(rhs.index(2));
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for Cartessian<T, 3>
            where T : Float + Scalar{
            fn $mth(&mut self, rhs : &'a Spherical<T>){
                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
                self[2].$mth(rhs.index(2));
            }
        }

        impl<T> $trt<CartessianND<T>> for Spherical<T>
            where T : Float + Scalar + AbsDiffEq<Epsilon = T>{
            fn $mth(&mut self, rhs : CartessianND<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);
                let z = self.index(2) $operate *rhs.index(2);

                self.radius = (x * x + y * y + z * z).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(z / self.radius).acos()};
                self.phi = y.atan2(x);
            }
        }

        impl<'a, T> $trt<&'a CartessianND<T>> for Spherical<T>
            where T : Float + Scalar + AbsDiffEq<Epsilon = T>{
            fn $mth(&mut self, rhs : &'a CartessianND<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                let x = self.index(0) $operate *rhs.index(0);
                let y = self.index(1) $operate *rhs.index(1);
                let z = self.index(2) $operate *rhs.index(2);

                self.radius = (x * x + y * y + z * z).sqrt();
                self.theta = if self.radius.abs_diff_eq(&T::zero(), <T as AbsDiffEq>::default_epsilon()) {T::zero()} else {(z / self.radius).acos()};
                self.phi = y.atan2(x);
            }
        }

        impl<T> $trt<Spherical<T>> for CartessianND<T>
            where T : Float + Scalar{
            fn $mth(&mut self, rhs : Spherical<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
                self[2].$mth(rhs.index(2));
            }
        }

        impl<'a, T> $trt<&'a Spherical<T>> for CartessianND<T>
            where T : Float + Scalar{
            fn $mth(&mut self, rhs : &'a Spherical<T>){
                if self.dim() != rhs.dim() {
                    panic!("Dimension of CartessianND should be 3 for operate with Spherical");
                }

                self[0].$mth(rhs.index(0));
                self[1].$mth(rhs.index(1));
                self[2].$mth(rhs.index(2));
            }
        }
    }
}

impl_assign_op!(AddAssign, +, add_assign, "");
impl_assign_op!(SubAssign, -, sub_assign, "");

impl<T> Neg for Spherical<T>
where
    T: FloatConst + Neg<Output = T> + Scalar + Rem<Output = T> + PartialOrd,
{
    type Output = Self;

    fn neg(mut self) -> Self {
        let pi = <T as FloatConst>::PI();
        self.theta = theta_rep(pi - self.theta);
        self.phi = (self.phi + pi) % (pi + pi);
        self
    }
}

impl<T> Dot<Spherical<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: Spherical<T>) -> T {
        self.radius
            * rhs.radius
            * (self.theta.sin() * rhs.theta.sin() * (self.phi - rhs.phi).cos()
                + self.theta.cos() * rhs.theta.cos())
    }
}

impl<'a, T> Dot<&'a Spherical<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: &'a Spherical<T>) -> T {
        self.radius
            * rhs.radius
            * (self.theta.sin() * rhs.theta.sin() * (self.phi - rhs.phi).cos()
                + self.theta.cos() * rhs.theta.cos())
    }
}

impl<T> Dot<Cartessian<T, 3>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: Cartessian<T, 3>) -> T {
        self.index(0) * *rhs.index(0)
            + self.index(1) * *rhs.index(1)
            + self.index(2) * *rhs.index(2)
    }
}

impl<'a, T> Dot<&'a Cartessian<T, 3>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: &'a Cartessian<T, 3>) -> T {
        self.index(0) * *rhs.index(0)
            + self.index(1) * *rhs.index(1)
            + self.index(2) * *rhs.index(2)
    }
}

impl<T> Dot<Spherical<T>> for Cartessian<T, 3>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: Spherical<T>) -> T {
        *self.index(0) * rhs.index(0)
            + *self.index(1) * rhs.index(1)
            + *self.index(2) * rhs.index(2)
    }
}

impl<'a, T> Dot<&'a Spherical<T>> for Cartessian<T, 3>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: &'a Spherical<T>) -> T {
        *self.index(0) * rhs.index(0)
            + *self.index(1) * rhs.index(1)
            + *self.index(2) * rhs.index(2)
    }
}

impl<T> Dot<CartessianND<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: CartessianND<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 3 for operate with Spherical");
        }

        self.index(0) * *rhs.index(0)
            + self.index(1) * *rhs.index(1)
            + self.index(2) * *rhs.index(2)
    }
}

impl<'a, T> Dot<&'a CartessianND<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: &'a CartessianND<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 3 for operate with Spherical");
        }

        self.index(0) * *rhs.index(0)
            + self.index(1) * *rhs.index(1)
            + self.index(2) * *rhs.index(2)
    }
}

impl<T> Dot<Spherical<T>> for CartessianND<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: Spherical<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 3 for operate with Spherical");
        }

        *self.index(0) * rhs.index(0)
            + *self.index(1) * rhs.index(1)
            + *self.index(2) * rhs.index(2)
    }
}

impl<'a, T> Dot<&'a Spherical<T>> for CartessianND<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn dot(&self, rhs: &'a Spherical<T>) -> T {
        if self.dim() != rhs.dim() {
            panic!("Dimension of CartessianND should be 3 for operate with Spherical");
        }

        *self.index(0) * rhs.index(0)
            + *self.index(1) * rhs.index(1)
            + *self.index(2) * rhs.index(2)
    }
}

impl<T> InnerProduct<Spherical<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: Spherical<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Spherical<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Spherical<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<Cartessian<T, 3>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: Cartessian<T, 3>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Cartessian<T, 3>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Cartessian<T, 3>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<Spherical<T>> for Cartessian<T, 3>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: Spherical<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Spherical<T>> for Cartessian<T, 3>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Spherical<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<CartessianND<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: CartessianND<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a CartessianND<T>> for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a CartessianND<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> InnerProduct<Spherical<T>> for CartessianND<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: Spherical<T>) -> T {
        self.dot(rhs)
    }
}

impl<'a, T> InnerProduct<&'a Spherical<T>> for CartessianND<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn inner_product(&self, rhs: &'a Spherical<T>) -> T {
        self.dot(rhs)
    }
}

impl<T> Norm for Spherical<T>
where
    T: Float + Scalar,
{
    type Output = T;

    fn norm_l1(&self) -> Self::Output {
        self.radius
            * (self.theta.cos().abs()
                + self.theta.sin().abs() * (self.phi.cos().abs() + self.phi.sin().abs()))
    }

    fn norm_l2(&self) -> Self::Output {
        self.radius
    }

    fn norm_l2_sqr(&self) -> Self::Output {
        self.radius
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
        let a = Spherical::<f64>::new(1.0f64, PI / 4f64, PI / 4f64);
        let b = Spherical::<f64>::new(1.0f64, PI * 1.75f64, PI / 4f64);
        let c = Spherical::<f64>::new(-1.0f64, PI * 1.25f64, PI * 1.25f64);

        assert_eq!(&a, &b);
        assert_eq!(&a, &c);
        a.abs_diff_eq(&b, 1e-6);
        a.abs_diff_eq(&c, 1e-6);
    }

    #[test]
    fn test_basic() {
        let mut a = Spherical::<f64>::new(1.0f64, PI / 4f64, PI / 4f64);
        assert_eq!(a.get_radius(), &1.0f64);
        assert_eq!(a.get_phi(), &(PI / 4f64));
        assert_eq!(a.get_theta(), &(PI / 4f64));

        *a.get_mut_radius() = 2.0f64;
        assert_abs_diff_eq!(a.radius, 2.0f64);

        *a.get_mut_phi() = PI / 2f64;
        assert_abs_diff_eq!(a.phi, PI / 2f64);

        *a.get_mut_theta() = PI / 2f64;
        assert_abs_diff_eq!(a.theta, PI / 2f64);

        assert_eq!(a.dim(), 3);

        let default = Spherical::<f64>::default();
        assert_eq!(
            default,
            Spherical::<f64> {
                radius: 0f64,
                theta: 0f64,
                phi: 0f64
            }
        );

        assert_abs_diff_eq!(a.index(0), 0.0f64);
        assert_abs_diff_eq!(a.index(1), 2.0f64);
        assert_abs_diff_eq!(a.index(2), 0.0f64);
    }

    #[test]
    fn test_from() {
        let a = Spherical::<f64>::new(1.0f64, PI / 4f64, PI / 4f64);
        assert_abs_diff_eq!(
            Cartessian::from(&a),
            Cartessian::new([0.5f64, 0.5f64, 0.5f64.sqrt()])
        );
        assert_abs_diff_eq!(
            CartessianND::from(&a),
            CartessianND::new(vec![0.5f64, 0.5f64, 0.5f64.sqrt()])
        );

        let carte = Cartessian::new([0.0f64, 1.0f64, 1f64]);
        assert_abs_diff_eq!(
            Spherical::<f64>::from(&carte),
            Spherical::<f64>::new(2.0f64.sqrt(), PI / 4f64, PI / 2f64)
        );

        let carte = CartessianND::new(vec![1.0f64, 0.0f64, 1f64]);
        assert_abs_diff_eq!(
            Spherical::<f64>::from(&carte),
            Spherical::<f64>::new(2.0f64.sqrt(), PI / 4f64, 0.0f64)
        );
    }

    #[test]
    #[should_panic]
    fn test_from_panic() {
        let carte = CartessianND::new(vec![1.0f64, 0.0f64]);
        let _x = Spherical::<f64>::from(&carte);
    }

    #[test]
    fn test_binary_op() {
        let mut a = Spherical::<f64>::new(1.0f64, 0f64, 0f64);
        let b = Spherical::<f64>::new(1.0f64, PI / 2f64, 0f64);
        let mut carte = Cartessian::new([-1f64, 0f64, 0f64]);
        let mut carte_nd = CartessianND::new(vec![-1f64, 0f64, 0f64]);

        (&a + &b).abs_diff_eq(
            &Spherical::<f64>::new(2.0f64.sqrt(), PI * 0.25f64, 0f64),
            1e-6,
        );
        (&a - &b).abs_diff_eq(
            &Spherical::<f64>::new(2.0f64.sqrt(), PI * 0.25f64, PI),
            1e-6,
        );
        (&b + &a).abs_diff_eq(
            &Spherical::<f64>::new(2.0f64.sqrt(), PI * 0.25f64, 0f64),
            1e-6,
        );
        (&b - &a).abs_diff_eq(
            &Spherical::<f64>::new(2.0f64.sqrt(), PI * 0.75f64, 0f64),
            1e-6,
        );

        (&a + &carte).abs_diff_eq(&Cartessian::new([-1f64, 0f64, 1f64]), 1e-6);
        (&a - &carte).abs_diff_eq(&Cartessian::new([1f64, 0f64, 1f64]), 1e-6);
        (&carte + &a).abs_diff_eq(&Cartessian::new([-1f64, 0f64, 1f64]), 1e-6);
        (&carte - &a).abs_diff_eq(&Cartessian::new([-1f64, 0f64, -1f64]), 1e-6);

        (&a + &carte_nd).abs_diff_eq(&CartessianND::new(vec![-1f64, 0f64, 1f64]), 1e-6);
        (&a - &carte_nd).abs_diff_eq(&CartessianND::new(vec![1f64, 0f64, 1f64]), 1e-6);
        (&carte_nd + &a).abs_diff_eq(&CartessianND::new(vec![-1f64, 0f64, 1f64]), 1e-6);
        (&carte_nd - &a).abs_diff_eq(&CartessianND::new(vec![-1f64, 0f64, -1f64]), 1e-6);

        (&a * 2.0f64).abs_diff_eq(&Spherical::<f64>::new(2f64, 0f64, 0f64), 1e-6);
        (2f64 * &a).abs_diff_eq(&Spherical::<f64>::new(2f64, 0f64, 0f64), 1e-6);
        (&a / 2.0f64).abs_diff_eq(&Spherical::<f64>::new(0.5f64, 0f64, 0f64), 1e-6);

        a += &b;
        a.abs_diff_eq(
            &Spherical::<f64>::new(2.0f64.sqrt(), PI * 0.25f64, 0f64),
            1e-6,
        );
        a -= &b;
        a.abs_diff_eq(&Spherical::<f64>::new(1.0f64, 0f64, 0f64), 1e-6);

        a += &carte;
        a.abs_diff_eq(
            &Spherical::<f64>::new(2.0f64.sqrt(), PI * 0.25f64, PI),
            1e-6,
        );
        a -= &carte;
        a.abs_diff_eq(&Spherical::<f64>::new(1.0f64, 0f64, 0f64), 1e-6);

        a += &carte_nd;
        a.abs_diff_eq(
            &Spherical::<f64>::new(2.0f64.sqrt(), PI * 0.25f64, PI),
            1e-6,
        );
        a -= &carte_nd;
        a.abs_diff_eq(&Spherical::<f64>::new(1.0f64, 0f64, 0f64), 1e-6);

        carte += a;
        carte.abs_diff_eq(&Cartessian::new([-1f64, 0f64, 1f64]), 1e-6);
        carte -= a;
        carte.abs_diff_eq(&Cartessian::new([-1f64, 0f64, 0f64]), 1e-6);

        carte_nd += a;
        carte_nd.abs_diff_eq(&CartessianND::new(vec![-1f64, 0f64, 1f64]), 1e-6);
        carte_nd -= a;
        carte_nd.abs_diff_eq(&CartessianND::new(vec![-1f64, 0f64, 0f64]), 1e-6);

        a = -a;
        a.abs_diff_eq(&Spherical::<f64>::new(1f64, PI, 0f64), 1e-6);
    }

    #[test]
    #[should_panic]
    fn test_add_panic() {
        let a = Spherical::<f64>::new(1.0f64, PI / 2f64, 0f64);
        let carte = CartessianND::new(vec![-1f64, 0f64]);

        let _x = carte + a;
    }

    #[test]
    fn test_dot() {
        let a = Spherical::<f64>::new(1.0f64, PI / 4f64, 0f64);
        let b = Spherical::<f64>::new(1.0f64, PI, 0f64);
        let carte = Cartessian::new([-1f64, 0f64, 0f64]);
        let carte_nd = CartessianND::new(vec![-1f64, 0f64, 0f64]);
        println!("{:?}", Cartessian::from(&a));
        println!("{:?}", Cartessian::from(&b));
        assert_abs_diff_eq!(a.dot(&b), -0.5f64.sqrt());
        assert_abs_diff_eq!(a.dot(&carte), -0.5f64.sqrt());
        assert_abs_diff_eq!(a.dot(&carte_nd), -0.5f64.sqrt());
        assert_abs_diff_eq!(carte.dot(&a), -0.5f64.sqrt());
        assert_abs_diff_eq!(carte_nd.dot(&a), -0.5f64.sqrt());
    }

    #[test]
    #[should_panic]
    fn test_dot_panic() {
        let a = Spherical::<f64>::new(1.0f64, PI / 2f64, 0f64);
        let carte = CartessianND::new(vec![-1f64, 0f64]);

        a.dot(&carte);
    }

    #[test]
    fn test_norm() {
        let a = Spherical::<f64>::new(1.0f64, PI / 4f64, 0f64);

        assert_abs_diff_eq!(a.norm_l1(), 2f64.sqrt());
        assert_abs_diff_eq!(a.norm_l2(), 1f64);
    }

    #[test]
    fn test_serde_spherical() {
        let a: Spherical<f64> = Spherical::new(2f64, 3.1415f64, 0f64);
        let expected = r#"{"radius":2.0,"theta":3.1415,"phi":0.0}"#;
        assert_eq!(expected, to_string(&a).unwrap());

        let expected: Spherical<f64> = from_str(&expected).unwrap();
        assert_eq!(a, expected);
    }
}
