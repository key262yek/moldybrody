
use num_traits::Float;
use num_traits::FloatConst;
use approx::AbsDiffEq;

use crate::vector::Vector;
use crate::{vector::arithmetic::Scalar, boundary::{plane::{SimplePlane, Plane, SimplePlanePair, PlanePair, SimpleBox, Parallelogram, Cube}, sphere::Sphere}};
use std::str::FromStr;
use crate::{error::{Error, ErrorCode}, vector::{Cartessian, CartessianND, curvilinear::{polar::Polar, spherical::Spherical}}};

pub trait FromArg : FromStr<Err = Error>{
    fn from_arg<I : Iterator<Item = String>>(iterable : &mut I) -> Result<Self, Error>{
        match iterable.next(){
            None => {
                return Err(Error::make_error_syntax(ErrorCode::InvalidNumberOfArguments));
            },
            Some(x) => {
                return x.parse::<Self>();
            }
        }
    }
}

impl<T : Scalar + Default + FromStr, const N : usize> FromArg for Cartessian<T, N> {}
impl<T : Scalar + Default + FromStr> FromArg for CartessianND<T> {}
impl<T : FromStr + Default + AbsDiffEq<Epsilon = T> + Float + FloatConst> FromArg for Polar<T> {}
impl<T : FromStr + Default + AbsDiffEq<Epsilon = T> + Float + FloatConst> FromArg for Spherical<T> {}

pub struct Brief<T>(pub T);

impl<T : Clone> Brief<T>{
    pub fn new(arg : &T) -> Self{
        Brief(arg.clone())
    }
}

pub trait ConvertBrief : Clone{
    fn brief(&self) -> Brief<Self>{
        Brief::new(self)
    }
}

impl<T : Scalar + Clone, const N : usize> ConvertBrief for Cartessian<T, N> {}
impl<T : Scalar + Clone> ConvertBrief for CartessianND<T> {}
impl<T : Clone> ConvertBrief for Polar<T> {}
impl<T : Clone> ConvertBrief for Spherical<T> {}
impl<T : Clone> ConvertBrief for SimplePlane<T> {}
impl<T : Clone> ConvertBrief for SimplePlanePair<T> {}
impl<V : Vector + Clone> ConvertBrief for Plane<V> {}
impl<V : Vector + Clone> ConvertBrief for PlanePair<V> {}
impl<T : Clone, const N : usize> ConvertBrief for SimpleBox<T, N> {}
impl<V : Vector + Clone> ConvertBrief for Cube<V> {}
impl<T : Scalar + Clone, const N : usize> ConvertBrief for Parallelogram<T, N> {}
impl<V : Vector + Clone> ConvertBrief for Sphere<V> {}

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::{Cartessian2D, Cartessian3D, Cartessian4D, CartessianND};

    #[test]
    fn test_from_arg(){
        let mut args = ["(2,3)", "(3,4)", "(4,5,6)", "(7,8,9,10)"].iter().map(|x| x.to_string());

        assert_eq!(Cartessian2D::<i32>::from_arg(&mut args).unwrap(), Cartessian2D::new([2, 3]));
        assert_eq!(Cartessian2D::<usize>::from_arg(&mut args).unwrap(), Cartessian2D::new([3usize, 4usize]));
        assert_eq!(Cartessian3D::<f32>::from_arg(&mut args).unwrap(), Cartessian3D::new([4f32, 5f32, 6f32]));
        assert_eq!(Cartessian4D::<f64>::from_arg(&mut args).unwrap(), Cartessian4D::new([7f64, 8f64, 9f64, 10f64]));
        assert_eq!(Cartessian4D::<f64>::from_arg(&mut args), Err(Error::make_error_syntax(ErrorCode::InvalidNumberOfArguments)));
    }

    #[test]
    fn test_from_arg_nd(){
        let mut args = ["(2,3)", "(3,4)", "(4,5,6)", "(7,8,9,10)"].iter().map(|x| x.to_string());

        assert_eq!(CartessianND::<i32>::from_arg(&mut args).unwrap(), CartessianND::new(vec![2, 3]));
        assert_eq!(CartessianND::<usize>::from_arg(&mut args).unwrap(), CartessianND::new(vec![3usize, 4usize]));
        assert_eq!(CartessianND::<f32>::from_arg(&mut args).unwrap(), CartessianND::new(vec![4f32, 5f32, 6f32]));
        assert_eq!(CartessianND::<f64>::from_arg(&mut args).unwrap(), CartessianND::new(vec![7f64, 8f64, 9f64, 10f64]));
        assert_eq!(CartessianND::<f64>::from_arg(&mut args), Err(Error::make_error_syntax(ErrorCode::InvalidNumberOfArguments)));
    }
}
