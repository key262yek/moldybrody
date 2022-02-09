


use crate::format_convert::{Brief, ConvertBrief};
use std::fmt::Formatter;
use std::fmt::{self, Display};
use std::fmt::LowerExp;
use std::ops::Rem;
use crate::{prelude::Error, vector::Dim};
use crate::boundary::Boundary;
use std::convert::TryInto;
use std::ops::Div;
use crate::vector::product::Norm;
use approx::AbsDiffEq;

use crate::vector::product::Dot;
use crate::vector::Vector;
use std::ops::Neg;
use crate::vector::arithmetic::Scalar;
use crate::vector::{Cartessian, CartessianND};

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Direction{
    Positive,
    Negative,
}

impl Display for Direction{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self{
            &Direction::Positive => {
                write!(f, "Positive")
            },
            &Direction::Negative => {
                write!(f, "Negative")
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SimplePlane<T>{
    idx : usize,
    pos : T,
    dir_in : Direction,
}

impl<T> SimplePlane<T>{
    pub fn new(idx : usize, pos : T, dir_in : Direction) -> Self{
        Self{
            idx,
            pos,
            dir_in,
        }
    }
}

macro_rules! impl_float_simpleplane {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for SimplePlane<$ty>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                match self.dir_in{
                    Direction::Positive => {
                        pos[self.idx] >= self.pos
                    },
                    Direction::Negative => {
                        pos[self.idx] <= self.pos
                    },
                }
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if self.pos.abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    let mut normal_vec = Cartessian::<$ty, N>::zeros();

                    match self.dir_in{
                        Direction::Positive => {
                            normal_vec[self.idx] = 1 as $ty;
                        },
                        Direction::Negative => {
                            normal_vec[self.idx] = -1 as $ty;
                        },
                    }

                    return Some(normal_vec);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, _pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N>{
                let mut normal_vec = Cartessian::<$ty, N>::zeros();

                match self.dir_in{
                    Direction::Positive => {
                        normal_vec[self.idx] = 1 as $ty;
                    },
                    Direction::Negative => {
                        normal_vec[self.idx] = -1 as $ty;
                    },
                }

                return normal_vec;
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos + movement;
                if self.check_inclusion(&destination){
                    return None;
                }

                match self.dir_in{
                    Direction::Positive => {
                        let t = (self.pos - pos[self.idx]) / movement[self.idx];
                        return Some(pos + movement * t);
                    },
                    Direction::Negative => {
                        let t = (self.pos - pos[self.idx]) / movement[self.idx];
                        return Some(pos + movement * t);
                    },
                }
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let destination = pos + movement;
                if self.check_inclusion(&destination){
                    return None;
                }

                match self.dir_in{
                    Direction::Positive => {
                        let t = (self.pos - pos[self.idx]) / movement[self.idx];
                        return Some(pos + movement * t);
                    },
                    Direction::Negative => {
                        let t = (self.pos - pos[self.idx]) / movement[self.idx];
                        return Some(pos + movement * t);
                    },
                }
            }
        }

        impl Boundary<CartessianND<$ty>> for SimplePlane<$ty>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }


                match self.dir_in{
                    Direction::Positive => {
                        pos[self.idx] >= self.pos
                    },
                    Direction::Negative => {
                        pos[self.idx] <= self.pos
                    },
                }
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if self.pos.abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    let mut normal_vec = CartessianND::zeros(pos.dim());

                    match self.dir_in{
                        Direction::Positive => {
                            normal_vec[self.idx] = 1 as $ty;
                        },
                        Direction::Negative => {
                            normal_vec[self.idx] = -1 as $ty;
                        },
                    }
                    return Some(normal_vec)
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty>{
                let mut normal_vec = CartessianND::zeros(pos.dim());

                match self.dir_in{
                    Direction::Positive => {
                        normal_vec[self.idx] = 1 as $ty;
                    },
                    Direction::Negative => {
                        normal_vec[self.idx] = -1 as $ty;
                    },
                }
                return normal_vec
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                match self.dir_in{
                    Direction::Positive => {
                        if pos[self.idx] + movement[self.idx] >= self.pos{
                            return None;
                        } else {
                            let t = (self.pos - pos[self.idx]) / movement[self.idx];
                            let mut result = pos.clone();
                            result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                            return Some(result);
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            let t = (self.pos - pos[self.idx]) / movement[self.idx];
                            let mut result = pos.clone();
                            result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                            return Some(result);
                        }
                    },
                }
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {

                match self.dir_in{
                    Direction::Positive => {
                        if pos[self.idx] + movement[self.idx] >= self.pos{
                            return None;
                        } else {
                            let t = (self.pos - pos[self.idx]) / movement[self.idx];
                            let mut result = pos.clone();
                            result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                            return Some(result);
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            let t = (self.pos - pos[self.idx]) / movement[self.idx];
                            let mut result = pos.clone();
                            result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                            return Some(result);
                        }
                    },
                }
            }
        }
    };
}

impl_float_simpleplane!(f32);
impl_float_simpleplane!(f64);

macro_rules! impl_int_simpleplane {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for SimplePlane<$ty>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                match self.dir_in{
                    Direction::Positive => {
                        pos[self.idx] >= self.pos
                    },
                    Direction::Negative => {
                        pos[self.idx] <= self.pos
                    },
                }
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if self.pos == pos[self.idx]{
                    let mut normal_vec = Cartessian::<$ty, N>::zeros();

                    match self.dir_in{
                        Direction::Positive => {
                            normal_vec[self.idx] = 1 as $ty;
                        },
                        Direction::Negative => {
                            normal_vec[self.idx] = -1 as $ty;
                        },
                    }

                    return Some(normal_vec);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, _pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N>{
                let mut normal_vec = Cartessian::<$ty, N>::zeros();

                match self.dir_in{
                    Direction::Positive => {
                        normal_vec[self.idx] = 1 as $ty;
                    },
                    Direction::Negative => {
                        normal_vec[self.idx] = -1 as $ty;
                    },
                }

                return normal_vec;
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                match self.dir_in{
                    Direction::Positive => {
                        if pos[self.idx] + movement[self.idx] >= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                }
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                match self.dir_in{
                    Direction::Positive => {
                        if pos[self.idx] + movement[self.idx] >= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                }
            }
        }

        impl Boundary<CartessianND<$ty>> for SimplePlane<$ty>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }


                match self.dir_in{
                    Direction::Positive => {
                        pos[self.idx] - self.pos >= 0 as $ty
                    },
                    Direction::Negative => {
                        pos[self.idx] - self.pos <= 0 as $ty
                    },
                }
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if self.pos == pos[self.idx]{
                    let mut normal_vec = CartessianND::zeros(pos.dim());

                    match self.dir_in{
                        Direction::Positive => {
                            normal_vec[self.idx] = 1 as $ty;
                        },
                        Direction::Negative => {
                            normal_vec[self.idx] = -1 as $ty;
                        },
                    }
                    return Some(normal_vec)
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty>{
                let mut normal_vec = CartessianND::zeros(pos.dim());

                match self.dir_in{
                    Direction::Positive => {
                        normal_vec[self.idx] = 1 as $ty;
                    },
                    Direction::Negative => {
                        normal_vec[self.idx] = -1 as $ty;
                    },
                }
                return normal_vec;
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                match self.dir_in{
                    Direction::Positive => {
                        if pos[self.idx] + movement[self.idx] >= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                }
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                match self.dir_in{
                    Direction::Positive => {
                        if pos[self.idx] + movement[self.idx] >= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            return Some(pos.clone());
                        }
                    },
                }
            }
        }
    };
}

impl_int_simpleplane!(i8);
impl_int_simpleplane!(i16);
impl_int_simpleplane!(i32);
impl_int_simpleplane!(i64);
impl_int_simpleplane!(i128);
impl_int_simpleplane!(isize);






impl<T : Display> Display for SimplePlane<T>{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.dir_in{
            Direction::Positive => write!(f, "SimplePlane at x[{}] = {} with positive to inside", self.idx, self.pos),
            Direction::Negative => write!(f, "SimplePlane at x[{}] = {} with negative to inside", self.idx, self.pos),
        }
    }
}

impl<T : LowerExp> LowerExp for SimplePlane<T>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result{
        match self.dir_in{
            Direction::Positive => {
                write!(f, "SimplePlane at x[{}] = ", self.idx)?;
                LowerExp::fmt(&self.pos, f)?;
                write!(f, " with positive to inside")
            },
            Direction::Negative => {
                write!(f, "SimplePlane at x[{}] = ", self.idx)?;
                LowerExp::fmt(&self.pos, f)?;
                write!(f, " with negative to inside")
            },
        }
    }
}



impl<T : Display> Display for Brief<SimplePlane<T>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.0.idx, self.0.pos, self.0.dir_in)
    }
}

impl<T : LowerExp> LowerExp for Brief<SimplePlane<T>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{},", self.0.idx)?;
        LowerExp::fmt(&self.0.pos, f)?;
        write!(f, ",{}", self.0.dir_in)
    }
}


#[derive(Clone, Debug, PartialEq)]
pub struct Plane<V : Vector>{
    normal_vec : V,  // direction : inside
    constant : <V as Vector>::Item,
}

impl<'a, V> Plane<V>
    where V : Vector + Clone + Norm + Div<<V as Norm>::Output, Output = V> + Dot<&'a V, Output = <V as Vector>::Item> + 'a,
          <V as Vector>::Item : Clone + Div<<V as Norm>::Output, Output = <V as Vector>::Item>,
          <V as Norm>::Output : Copy{
    pub fn new(normal_vec : &V, constant : <V as Vector>::Item) -> Self{
        let n = normal_vec.norm_l2();
        Self{
            normal_vec : (*normal_vec).clone() / n,
            constant : constant / n,
        }
    }

    pub fn from_point(normal_vec : &'a V, point : &'a V) -> Self{
        let n = normal_vec.norm_l2();
        let vec = normal_vec.clone() / n;
        let c = vec.dot(point);
        Self{
            normal_vec : vec,
            constant : c / n,
        }
    }
}

macro_rules! impl_float_plane {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for Plane<Cartessian<$ty, N>>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                let d = self.normal_vec.dot(pos);

                d >= self.constant
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                let d = self.normal_vec.dot(pos);
                if self.constant.abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return Some(self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, _pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N>{
                self.normal_vec.clone()
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos + movement;
                if self.check_inclusion(&destination){
                    return None;
                }

                let t = - (self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(pos + movement * t);
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                let destination = pos + movement;
                if self.check_inclusion(&destination){
                    return None;
                }

                let t = - (self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(pos + movement * t);
            }
        }

        impl Boundary<CartessianND<$ty>> for Plane<CartessianND<$ty>>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                let d = self.normal_vec.dot(pos);

                d >= self.constant
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                let d = self.normal_vec.dot(pos);
                if self.constant.abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return Some(self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, _pos : &CartessianND<$ty>) -> CartessianND<$ty>{
                self.normal_vec.clone()
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{

                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let mut result = pos.clone();
                result.zip_mut_with(movement, move |x, y| *x = *x + *y);
                if self.check_inclusion(&result){
                    return None;
                }

                let t = - (self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                result.clone_from(pos);
                result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                return Some(result);
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                let mut result = pos.clone();
                result.zip_mut_with(movement, move |x, y| *x = *x + *y);
                if self.check_inclusion(&result){
                    return None;
                }

                let t = - (self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                result.clone_from(pos);
                result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                return Some(result);
            }
        }
    };
}

impl_float_plane!(f32);
impl_float_plane!(f64);


impl<T, V> Display for Plane<V>
    where V : Vector<Item = T> + Display + ConvertBrief,
          Brief<V> : Display,
          T : Display{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Plane normal to ({}) with constant {}", self.normal_vec.brief(), self.constant)
    }
}

impl<T, V> LowerExp for Plane<V>
    where V : Vector<Item = T> + LowerExp + ConvertBrief,
          Brief<V> : LowerExp,
          T : LowerExp{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result{
        write!(f, "Plane normal to (")?;
        LowerExp::fmt(&self.normal_vec.brief(), f)?;
        write!(f, ") with constant ")?;
        LowerExp::fmt(&self.constant, f)
    }
}

impl<V, T> Display for Brief<Plane<V>>
    where V : Vector<Item = T> + Display + ConvertBrief,
          Brief<V> : Display,
          T : Display{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{},{}", &self.0.normal_vec.brief(), &self.0.constant)
    }
}

impl<V, T> LowerExp for Brief<Plane<V>>
    where V : Vector<Item = T> + LowerExp + ConvertBrief,
          Brief<V> : LowerExp,
          T : LowerExp{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        LowerExp::fmt(&self.0.normal_vec.brief(), f)?;
        write!(f, ",")?;
        LowerExp::fmt(&self.0.constant, f)
    }
}


#[derive(Clone, Debug, PartialEq)]
pub struct SimplePlanePair<T>{
    idx : usize,
    pos : [T; 2],
}

impl<T : PartialOrd + Copy + AbsDiffEq> SimplePlanePair<T>{
    pub fn new(idx : usize, pos : [T; 2]) -> Option<Self>{
        if pos[0].abs_diff_eq(&pos[1], <T as AbsDiffEq>::default_epsilon()){
            return None;
        }

        let (a, b) = match pos[0] > pos[1]{
            true => (pos[1], pos[0]),
            false => (pos[0], pos[1]),
        };

        Some(Self{
            idx,
            pos : [a, b],
        })
    }
}

macro_rules! impl_float_simpleplanepair {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for SimplePlanePair<$ty>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                self.pos[0] <= pos[self.idx]  && pos[self.idx] <= self.pos[1]
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut normal_vec = Cartessian::<$ty, N>::zeros();
                if self.pos[0].abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    normal_vec[self.idx] = 1 as $ty;
                    return Some(normal_vec);
                } else if self.pos[1].abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    normal_vec[self.idx] = -1 as $ty;
                    return Some(normal_vec);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N>{
                let mut normal_vec = Cartessian::<$ty, N>::zeros();
                if self.pos[0].abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    normal_vec[self.idx] = 1 as $ty;
                    return normal_vec;
                } else {
                    normal_vec[self.idx] = -1 as $ty;
                    return normal_vec;
                }
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    return Some(pos + movement * t);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    return Some(pos + movement * t);
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    return Some(pos + movement * t);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    return Some(pos + movement * t);
                } else {
                    return None;
                }
            }
        }

        impl Boundary<CartessianND<$ty>> for SimplePlanePair<$ty>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                self.pos[0] <= pos[self.idx]  && pos[self.idx] <= self.pos[1]
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut normal_vec = CartessianND::zeros(pos.dim());
                if self.pos[0].abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    normal_vec[self.idx] = 1 as $ty;
                    return Some(normal_vec);
                } else if self.pos[1].abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    normal_vec[self.idx] = -1 as $ty;
                    return Some(normal_vec);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty>{
                let mut normal_vec = CartessianND::zeros(pos.dim());
                if self.pos[0].abs_diff_eq(&pos[self.idx], <$ty as AbsDiffEq>::default_epsilon()){
                    normal_vec[self.idx] = 1 as $ty;
                    return normal_vec;
                } else {
                    normal_vec[self.idx] = -1 as $ty;
                    return normal_vec;
                }
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    let mut result = pos.clone();
                    result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                    return Some(result);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    let mut result = pos.clone();
                    result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                    return Some(result);
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    let mut result = pos.clone();
                    result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                    return Some(result);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    let mut result = pos.clone();
                    result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                    return Some(result);
                } else {
                    return None;
                }
            }
        }
    };
}

impl_float_simpleplanepair!(f32);
impl_float_simpleplanepair!(f64);

macro_rules! impl_int_simpleplanepair {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for SimplePlanePair<$ty>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                self.pos[0] <= pos[self.idx]  && pos[self.idx] <= self.pos[1]
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut normal_vec = Cartessian::<$ty, N>::zeros();
                if self.pos[0] == pos[self.idx]{
                    normal_vec[self.idx] = 1 as $ty;
                    return Some(normal_vec);
                } else if self.pos[1] == pos[self.idx]{
                    normal_vec[self.idx] = -1 as $ty;
                    return Some(normal_vec);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N>{
                let mut normal_vec = Cartessian::<$ty, N>::zeros();
                if self.pos[0] == pos[self.idx]{
                    normal_vec[self.idx] = 1 as $ty;
                    return normal_vec;
                } else {
                    normal_vec[self.idx] = -1 as $ty;
                    return normal_vec;
                }
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0] || self.pos[1] < destination{
                    return Some(pos.clone());
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0] || self.pos[1] < destination{
                    return Some(pos.clone());
                } else {
                    return None;
                }
            }
        }

        impl Boundary<CartessianND<$ty>> for SimplePlanePair<$ty>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                self.pos[0] <= pos[self.idx]  && pos[self.idx] <= self.pos[1]
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut normal_vec = CartessianND::zeros(pos.dim());
                if self.pos[0] == pos[self.idx]{
                    normal_vec[self.idx] = 1 as $ty;
                    return Some(normal_vec);
                } else if self.pos[1] == pos[self.idx]{
                    normal_vec[self.idx] = -1 as $ty;
                    return Some(normal_vec);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty>{
                let mut normal_vec = CartessianND::zeros(pos.dim());
                if self.pos[0] == pos[self.idx]{
                    normal_vec[self.idx] = 1 as $ty;
                    return normal_vec;
                } else {
                    normal_vec[self.idx] = -1 as $ty;
                    return normal_vec;
                }
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0] || self.pos[1] < destination{
                    return Some(pos.clone());
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0] || self.pos[1] < destination{
                    return Some(pos.clone());
                } else {
                    return None;
                }
            }
        }
    };
}

impl_int_simpleplanepair!(i8);
impl_int_simpleplanepair!(i16);
impl_int_simpleplanepair!(i32);
impl_int_simpleplanepair!(i64);
impl_int_simpleplanepair!(i128);
impl_int_simpleplanepair!(isize);


impl<T : Display> Display for SimplePlanePair<T>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "SimplePlanePair between x[{}] = {} and {}", self.idx, self.pos[0], self.pos[1])
    }
}

impl<T : LowerExp> LowerExp for SimplePlanePair<T>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result{
        write!(f, "SimplePlanePair between x[{}] = ", self.idx)?;
        LowerExp::fmt(&self.pos[0], f)?;
        write!(f, " and ")?;
        LowerExp::fmt(&self.pos[1], f)
    }
}


impl<T : Display> Display for Brief<SimplePlanePair<T>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", &self.0.idx, &self.0.pos[0], &self.0.pos[1])
    }
}

impl<T : LowerExp> LowerExp for Brief<SimplePlanePair<T>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{},", &self.0.idx)?;
        LowerExp::fmt(&self.0.pos[0], f)?;
        write!(f, ",")?;
        LowerExp::fmt(&self.0.pos[1], f)
    }
}




#[derive(Clone, Debug, PartialEq)]
pub struct PlanePair<V : Vector>{
    normal_vec : V,
    constant : [<V as Vector>::Item; 2],
}

impl<V> PlanePair<V>
    where V : Vector + Clone + Norm + Div<<V as Norm>::Output, Output = V>,
         <V as Vector>::Item : PartialOrd + AbsDiffEq + Div<<V as Norm>::Output, Output = <V as Vector>::Item>,
         <V as Norm>::Output : Copy{
    pub fn new(normal_vec : V, constant : [<V as Vector>::Item; 2]) -> Option<Self>{
        if constant[0].abs_diff_eq(&constant[1], <<V as Vector>::Item as AbsDiffEq>::default_epsilon()){
            return None;
        }

        let n = normal_vec.norm_l2();

        let (a, b) = match constant[0] > constant[1]{
            true => (constant[1], constant[0]),
            false => (constant[0], constant[1]),
        };
        Some(Self{
            normal_vec : normal_vec.clone() / n,
            constant : [a / n, b / n],
        })
    }
}

macro_rules! impl_float_planepair {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for PlanePair<Cartessian<$ty, N>>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                let d = self.normal_vec.dot(pos);

                self.constant[0] <= d && d <= self.constant[1]
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return Some(self.normal_vec.clone());
                } else if self.constant[1].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return Some(-self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N>{
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return self.normal_vec.clone();
                } else {
                    return -self.normal_vec.clone();
                }
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let d = self.normal_vec.dot(pos + movement);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                let d = self.normal_vec.dot(pos + movement);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else {
                    return None;
                }
            }
        }

        impl Boundary<CartessianND<$ty>> for PlanePair<CartessianND<$ty>>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                let d = self.normal_vec.dot(pos);

                self.constant[0] <= d && d <= self.constant[1]
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return Some(self.normal_vec.clone());
                } else if self.constant[1].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return Some(-self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty>{
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()){
                    return self.normal_vec.clone();
                } else {
                    return -self.normal_vec.clone();
                }
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let mut result = pos.clone();
                result.zip_mut_with(movement, |x, y| *x = *x + *y);
                let d = self.normal_vec.dot(&result);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                let mut result = pos.clone();
                result.zip_mut_with(movement, |x, y| *x = *x + *y);
                let d = self.normal_vec.dot(&result);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos)) / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else {
                    return None;
                }
            }
        }
    };
}

impl_float_planepair!(f32);
impl_float_planepair!(f64);

impl<T, V> Display for PlanePair<V>
    where V : Vector<Item = T> + Display + ConvertBrief,
          Brief<V> : Display,
          T : Display{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PlanePair normal to ({}) between constant {} and {}", self.normal_vec.brief(), self.constant[0], self.constant[1])
    }
}

impl<T, V> LowerExp for PlanePair<V>
    where V : Vector<Item = T> + LowerExp+ ConvertBrief,
          Brief<V> : LowerExp,
          T : LowerExp{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result{
        write!(f, "PlanePair normal to (")?;
        LowerExp::fmt(&self.normal_vec.brief(), f)?;
        write!(f, ") between constant ")?;
        LowerExp::fmt(&self.constant[0], f)?;
        write!(f, " and ")?;
        LowerExp::fmt(&self.constant[1], f)
    }
}


impl<V, T> Display for Brief<PlanePair<V>>
    where V : Vector<Item = T> + Display+ ConvertBrief,
          Brief<V> : Display,
          T : Display{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", &self.0.normal_vec.brief(), &self.0.constant[0], &self.0.constant[1])
    }
}

impl<V, T> LowerExp for Brief<PlanePair<V>>
    where V : Vector<Item = T> + LowerExp+ ConvertBrief,
          Brief<V> : LowerExp,
          T : LowerExp{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        LowerExp::fmt(&self.0.normal_vec.brief(), f)?;
        write!(f, ",")?;
        LowerExp::fmt(&self.0.constant[0], f)?;
        write!(f, ",")?;
        LowerExp::fmt(&self.0.constant[1], f)
    }
}


#[derive(Clone, Debug, PartialEq)]
pub struct SimpleBox<T, const N : usize>{
    planes : [SimplePlanePair<T>; N],
}

impl<T, const N : usize> SimpleBox<T, N>{
    pub fn new(planes : [SimplePlanePair<T>; N]) -> Result<Self, Error>
        where T : Clone{
        use itertools::Itertools;

        let count = planes.iter().map(|x| x.idx)
                               .unique().count();
        if count != N {
            return Err(Error::make_error_syntax(crate::prelude::ErrorCode::InvalidArgumentInput));
        }

        Ok(Self{
            planes : planes.clone()
        })
    }

    pub fn cube_with_center<'a, V>(center : &'a V, half : T) -> Result<Self, Error>
        where V : Vector<Item = T> + Dim<N>,
              &'a V : IntoIterator<Item = &'a T>,
              T : Scalar + PartialOrd + AbsDiffEq + Rem<Output = T>{

        if center.dim() != N  || half < T::zero(){
            return Err(Error::make_error_syntax(crate::prelude::ErrorCode::InvalidDimension));
        }

        let planes : [SimplePlanePair<T>; N]
                = (0..N).zip(center.into_iter())
                        .map(|(idx, y)| SimplePlanePair::new(idx, [*y - half, *y + half]).unwrap())
                        .collect::<Vec<SimplePlanePair<T>>>()
                        .try_into().unwrap();
        Ok(Self{
            planes,
        })
    }
}

impl<T, const N : usize> Boundary<Cartessian<T, N>> for SimpleBox<T, N>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
          SimplePlanePair<T> : Boundary<Cartessian<T, N>>{
    fn check_inclusion(&self, pos : &Cartessian<T, N>) -> bool{
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos : &Cartessian<T, N>) -> Option<Cartessian<T, N>>{
        for plane in &self.planes{
            match plane.normal_at(pos){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos : &Cartessian<T, N>) -> Cartessian<T, N> {
      for plane in &self.planes{
            match plane.normal_at(pos){
                Some(v) => {return v},
                None => {},
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(&self, pos : &Cartessian<T, N>, movement : &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        if !self.check_inclusion(pos){
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)){
            return None;
        }

        for plane in &self.planes{
            match plane.find_intersect_unsafe(pos, movement){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }

    fn find_intersect_unsafe(&self, pos : &Cartessian<T, N>, movement : &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        for plane in &self.planes{
            match plane.find_intersect_unsafe(pos, movement){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }
}

impl<T, const N : usize> Boundary<CartessianND<T>> for SimpleBox<T, N>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
          SimplePlanePair<T> : Boundary<CartessianND<T>>{
    fn check_inclusion(&self, pos : &CartessianND<T>) -> bool{
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos : &CartessianND<T>) -> Option<CartessianND<T>>{
        for plane in &self.planes{
            match plane.normal_at(pos){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos : &CartessianND<T>) -> CartessianND<T> {
      for plane in &self.planes{
            match plane.normal_at(pos){
                Some(v) => {return v},
                None => {},
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(&self, pos : &CartessianND<T>, movement : &CartessianND<T>) -> Option<CartessianND<T>> {
        if !self.check_inclusion(pos){
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)){
            return None;
        }

        for plane in &self.planes{
            match plane.find_intersect_unsafe(pos, movement){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }

    fn find_intersect_unsafe(&self, pos : &CartessianND<T>, movement : &CartessianND<T>) -> Option<CartessianND<T>> {
        for plane in &self.planes{
            match plane.find_intersect_unsafe(pos, movement){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }
}

impl<T : Clone + Display, const N : usize> Display for SimpleBox<T, N>{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SimpleBox enveloped in plane pairs ({}", &self.planes[0].brief())?;
        for plane in self.planes.iter().skip(1){
            write!(f, "), ({}", plane.brief())?;
        }
        write!(f, ")")
    }
}

impl<T : Clone + Display +  LowerExp, const N : usize> LowerExp for SimpleBox<T, N>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result{
        write!(f, "SimpleBox enveloped in plane pairs (")?;
        LowerExp::fmt(&self.planes[0].brief(), f)?;
        for plane in self.planes.iter().skip(1){
            write!(f, "), (")?;
            LowerExp::fmt(&plane.brief(), f)?;
        }
        write!(f, ")")
    }
}


impl<T : Clone + Display, const N : usize> Display for Brief<SimpleBox<T, N>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "({}", &self.0.planes[0].brief())?;
        for plane in self.0.planes.iter().skip(1){
            write!(f, "),({}", plane.brief())?;
        }
        write!(f, ")")
    }
}

impl<T : Clone + LowerExp, const N : usize> LowerExp for Brief<SimpleBox<T, N>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "(")?;
        LowerExp::fmt(&self.0.planes[0].brief(), f)?;
        for plane in self.0.planes.iter().skip(1){
            write!(f, "),(", )?;
            LowerExp::fmt(&plane.brief(), f)?;
        }
        write!(f, ")")
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Cube<V : Vector>{
    center : V,
    radius : <V as Vector>::Item,
}

impl<V : Vector> Cube<V>{
    pub fn new(center : &V, radius : <V as Vector>::Item) -> Self
        where V : Clone, <V as Vector>::Item : Clone{
        Self{
            center : center.clone(),
            radius : radius.clone(),
        }
    }
}

macro_rules! impl_float_cube {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for Cube<Cartessian<$ty, N>>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                for (c, x) in self.center.into_iter().zip(pos){
                    if (*x - *c) > self.radius || (*x - *c) < -self.radius{
                        return false
                    }
                }
                return true;
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                let mut vec = Cartessian::<$ty, N>::default();
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon){
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon){
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let mut vec = Cartessian::<$ty, N>::default();
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon){
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon){
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                for idx in 0..N{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius{
                        let t = (c - self.radius - p) / m;
                        return Some(pos + movement * t);
                    } else if p + m > c + self.radius{
                        let t = (c + self.radius - p) / m;
                        return Some(pos + movement * t);
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                for idx in 0..N{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius{
                        let t = (c - self.radius - p) / m;
                        return Some(pos + movement * t);
                    } else if p + m > c + self.radius{
                        let t = (c + self.radius - p) / m;
                        return Some(pos + movement * t);
                    }
                }
                return None;
            }
        }

        impl Boundary<CartessianND<$ty>> for Cube<CartessianND<$ty>>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                if self.center.dim() != pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                for (c, x) in self.center.into_iter().zip(pos){
                    if *x > *c + self.radius || *x < *c -self.radius{
                        return false
                    }
                }
                return true;
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                if self.center.dim() != pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon){
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon){
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty> {
                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon){
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon){
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let n = pos.dim();
                for idx in 0..n{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius{
                        let t = (c - self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(&movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    } else if p + m > c + self.radius{
                        let t = (c + self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(&movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                let n = pos.dim();
                for idx in 0..n{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius{
                        let t = (c - self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(&movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    } else if p + m > c + self.radius{
                        let t = (c + self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(&movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    }
                }
                return None;
            }
        }
    };
}

impl_float_cube!(f32);
impl_float_cube!(f64);

macro_rules! impl_int_cube {
    ($ty : ident) => {
        impl<const N : usize> Boundary<Cartessian<$ty, N>> for Cube<Cartessian<$ty, N>>{
            fn check_inclusion(&self, pos : &Cartessian<$ty, N>) -> bool{
                for (c, x) in self.center.into_iter().zip(pos){
                    if (*x - *c) > self.radius || (*x - *c) < -self.radius{
                        return false
                    }
                }
                return true;
            }

            fn normal_at(&self, pos : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>>{
                let mut vec = Cartessian::<$ty, N>::default();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if *x == *c + &self.radius{
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if *x == *c - &self.radius{
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let mut vec = Cartessian::<$ty, N>::default();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if *x == *c + &self.radius{
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if *x == *c - &self.radius{
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                for idx in 0..N{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius || p + m > c + self.radius{
                        return Some(pos.clone());
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                for idx in 0..N{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius || p + m > c + self.radius{
                        return Some(pos.clone());
                    }
                }

                return None;
            }
        }

        impl Boundary<CartessianND<$ty>> for Cube<CartessianND<$ty>>{
            fn check_inclusion(&self, pos : &CartessianND<$ty>) -> bool{
                if self.center.dim() != pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                for (c, x) in self.center.into_iter().zip(pos){
                    if *x > *c + self.radius || *x < *c -self.radius{
                        return false
                    }
                }
                return true;
            }

            fn normal_at(&self, pos : &CartessianND<$ty>) -> Option<CartessianND<$ty>>{
                if self.center.dim() != pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if *x == *c + &self.radius{
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if *x == *c - &self.radius{
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty> {
                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate(){
                    if *x == *c + &self.radius{
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if *x == *c - &self.radius{
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let n = pos.dim();
                for idx in 0..n{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius || p + m > c + self.radius{
                        return Some(pos.clone());
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                let n = pos.dim();
                for idx in 0..n{
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c -self.radius || p + m > c + self.radius{
                        return Some(pos.clone());
                    }
                }
                return None;
            }
        }
    };
}

impl_int_cube!(i8);
impl_int_cube!(i16);
impl_int_cube!(i32);
impl_int_cube!(i64);
impl_int_cube!(i128);
impl_int_cube!(isize);

impl<V, T> Display for Cube<V>
    where V : Vector<Item = T> + ConvertBrief + Display,
          Brief<V> : Display,
          T : Display{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "Cube has center at ({}) with radius {}", self.center.brief(), self.radius)
    }
}

impl<V, T> LowerExp for Cube<V>
    where V : Vector<Item = T> + ConvertBrief + LowerExp,
          Brief<V> : LowerExp,
          T : LowerExp{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result{
        write!(f, "Cube has center at (")?;
        LowerExp::fmt(&self.center.brief(), f)?;
        write!(f, ") with radius ")?;
        LowerExp::fmt(&self.radius, f)
    }
}

impl<V, T> Display for Brief<Cube<V>>
    where V : Vector<Item = T> + ConvertBrief + Display,
          Brief<V> : Display,
          T : Display{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{},{}", &self.0.center.brief(), &self.0.radius)
    }
}

impl<V, T> LowerExp for Brief<Cube<V>>
    where V : Vector<Item = T> + ConvertBrief + LowerExp,
          Brief<V> : LowerExp,
          T : LowerExp{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        LowerExp::fmt(&self.0.center.brief(), f)?;
        write!(f, ",")?;
        LowerExp::fmt(&self.0.radius, f)
    }
}


#[derive(Clone, Debug, PartialEq)]
pub struct Parallelogram<T : Scalar, const N : usize>{
    planes : [PlanePair<Cartessian<T, N>>; N],
}

impl<T : Scalar, const N : usize> Parallelogram<T, N>{
    pub fn new(planes : [PlanePair<Cartessian<T, N>>; N]) -> Result<Self, Error>
        where T : Clone + AbsDiffEq,
             <T as AbsDiffEq>::Epsilon : Clone{

        for i in 0..N-1 {
            for j in i + 1..N{
                if planes[i].normal_vec.abs_diff_eq(&planes[j].normal_vec, <T as AbsDiffEq>::default_epsilon()){
                    return Err(Error::make_error_syntax(crate::prelude::ErrorCode::InvalidArgumentInput));
                }
            }
        }

        Ok(Self{
            planes : planes.clone()
        })
    }
}

impl<T, const N : usize> Boundary<Cartessian<T, N>> for Parallelogram<T, N>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
          PlanePair<Cartessian<T, N>> : Boundary<Cartessian<T, N>>{
    fn check_inclusion(&self, pos : &Cartessian<T, N>) -> bool{
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos : &Cartessian<T, N>) -> Option<Cartessian<T, N>>{
        for plane in &self.planes{
            match plane.normal_at(pos){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos : &Cartessian<T, N>) -> Cartessian<T, N> {
        for plane in &self.planes{
            match plane.normal_at(pos){
                Some(v) => {return v},
                None => {},
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(&self, pos : &Cartessian<T, N>, movement : &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        if !self.check_inclusion(pos){
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)){
            return None;
        }

        for plane in &self.planes{
            match plane.find_intersect_unsafe(pos, movement){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }

    fn find_intersect_unsafe(&self, pos : &Cartessian<T, N>, movement : &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        for plane in &self.planes{
            match plane.find_intersect_unsafe(pos, movement){
                Some(v) => {return Some(v)},
                None => {},
            }
        }
        return None;
    }
}

impl<T : Display + LowerExp + Scalar, const N : usize> Display for Parallelogram<T, N>{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Parallelogram enveloped in plane pairs ({}", &self.planes[0].brief())?;
        for plane in self.planes.iter().skip(1){
            write!(f, "), ({}", plane.brief())?;
        }
        write!(f, ")")
    }
}

impl<T : Display + LowerExp + Scalar, const N : usize> LowerExp for Parallelogram<T, N>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result{
        write!(f, "Parallelogram enveloped in plane pairs (")?;
        LowerExp::fmt(&self.planes[0].brief(), f)?;
        for plane in self.planes.iter().skip(1){
            write!(f, "), (")?;
            LowerExp::fmt(&plane.brief(), f)?;
        }
        write!(f, ")")
    }
}


impl<T : Clone + Display + Scalar, const N : usize> Display for Brief<Parallelogram<T, N>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "({}", &self.0.planes[0].brief())?;
        for plane in self.0.planes.iter().skip(1){
            write!(f, "),({}", plane.brief())?;
        }
        write!(f, ")")
    }
}

impl<T : Clone + LowerExp + Scalar + Display, const N : usize> LowerExp for Brief<Parallelogram<T, N>>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "(")?;
        LowerExp::fmt(&self.0.planes[0].brief(), f)?;
        for plane in self.0.planes.iter().skip(1){
            write!(f, "),(", )?;
            LowerExp::fmt(&plane.brief(), f)?;
        }
        write!(f, ")")
    }
}



#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::{Cartessian2D, Cartessian3D, CartessianND};
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_simpleplane(){
        let plane = SimplePlane::new(0, 1f64, Direction::Positive);

        // float 2D
        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(Cartessian2D::new([1f64, 0.0])));
        assert_eq!(plane.normal_at_unsafe(&a), Cartessian2D::new([1f64, 0.0]));

        a[0] = 2f64;
        let mut movement = Cartessian2D::new([-2f64; 2]);
        assert_eq!(plane.find_intersect(&a, &movement), Some(Cartessian2D::new([1f64, 1f64])));
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([1f64, 1f64])));

        movement[0] = -0.5f64;
        assert_eq!(plane.find_intersect(&a, &movement), None);
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), None);

        // float ND
        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(CartessianND::new(vec![1f64, 0.0])));
        assert_eq!(plane.normal_at_unsafe(&a), CartessianND::new(vec![1f64, 0.0]));

        a[0] = 2f64;
        let mut movement = CartessianND::new(vec![-2f64; 2]);
        assert_eq!(plane.find_intersect(&a, &movement), Some(CartessianND::new(vec![1f64, 1f64])));
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), Some(CartessianND::new(vec![1f64, 1f64])));

        movement[0] = -0.5f64;
        assert_eq!(plane.find_intersect(&a, &movement), None);
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), None);

        let plane = SimplePlane::new(0, 1, Direction::Positive);
        // int 2D
        let mut a = Cartessian2D::new([2; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 0;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 1;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(Cartessian2D::new([1, 0])));
        assert_eq!(plane.normal_at_unsafe(&a), Cartessian2D::new([1, 0]));

        a[0] = 1;
        let mut movement = Cartessian2D::new([-1, 0]);
        assert_eq!(plane.find_intersect(&a, &movement), Some(Cartessian2D::new([1, 2])));
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([1, 2])));

        movement[0] = 0;
        assert_eq!(plane.find_intersect(&a, &movement), None);
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), None);

        // int ND
        let mut a = CartessianND::new(vec![2; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 0;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 1;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(CartessianND::new(vec![1, 0])));
        assert_eq!(plane.normal_at_unsafe(&a), CartessianND::new(vec![1, 0]));

        a[0] = 1;
        let mut movement = CartessianND::new(vec![-1, 0]);
        assert_eq!(plane.find_intersect(&a, &movement), Some(CartessianND::new(vec![1, 2])));
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), Some(CartessianND::new(vec![1, 2])));

        movement[0] = 0;
        assert_eq!(plane.find_intersect(&a, &movement), None);
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), None);


        // display
        assert_eq!(format!("{}", plane), "SimplePlane at x[0] = 1 with positive to inside");
        assert_eq!(format!("{:.3e}", plane), "SimplePlane at x[0] = 1.000e0 with positive to inside");
        assert_eq!(format!("{}", plane.brief()), "0,1,Positive");
        assert_eq!(format!("{:.3e}", plane.brief()), "0,1.000e0,Positive");
    }

    #[test]
    #[should_panic]
    fn test_simpleplane_panic(){
        let plane = SimplePlane::new(3, 1f64, Direction::Positive);

        let a = Cartessian2D::new([2f64; 2]);
        plane.check_inclusion(&a);
    }

    #[test]
    #[should_panic]
    fn test_simpleplane_find_itersect(){
        let plane = SimplePlane::new(0, 1f64, Direction::Positive);

        let a = Cartessian2D::new([0f64; 2]);
        let movement = Cartessian2D::new([1f64, 0f64]);
        plane.find_intersect(&a, &movement);
    }



    #[test]
    fn test_plane(){
        // 2d
        let plane = Plane::new(&Cartessian2D::new([1f64, 1f64]), 0f64);
        assert_abs_diff_eq!(plane.normal_vec, Cartessian2D::new([1.0 / 2f64.sqrt(); 2]));

        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -4f64;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -2f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(Cartessian2D::new([1.0 / 2f64.sqrt(); 2])));
        assert_eq!(plane.normal_at_unsafe(&a), Cartessian2D::new([1.0 / 2f64.sqrt(); 2]));

        a[0] = 2f64;
        let movement = Cartessian2D::new([-3f64; 2]);
        assert_eq!(plane.find_intersect(&a, &movement), Some(Cartessian2D::new([0f64; 2])));
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([0f64; 2])));


        // nd
        let plane = Plane::new(&CartessianND::new(vec![1f64, 1f64]), 0f64);
        assert_abs_diff_eq!(plane.normal_vec, CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]));

        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -4f64;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -2f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])));
        assert_eq!(plane.normal_at_unsafe(&a), CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]));

        a[0] = 2f64;
        let movement = CartessianND::new(vec![-3f64; 2]);
        assert_eq!(plane.find_intersect(&a, &movement), Some(CartessianND::new(vec![0f64; 2])));
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), Some(CartessianND::new(vec![0f64; 2])));

        assert_eq!(format!("{}", plane), "Plane normal to (0.7071067811865475:0.7071067811865475) with constant 0");
        assert_eq!(format!("{:.3e}", plane), "Plane normal to (7.071e-1:7.071e-1) with constant 0.000e0");
        assert_eq!(format!("{}", plane.brief()), "0.7071067811865475:0.7071067811865475,0");
        assert_eq!(format!("{:.3e}", plane.brief()), "7.071e-1:7.071e-1,0.000e0");
    }

    #[test]
    #[should_panic]
    fn test_plane_panic(){
        let plane = Plane::new(&CartessianND::new(vec![3f64; 3]), 0f64);

        let a = CartessianND::new(vec![2f64; 2]);
        plane.check_inclusion(&a);
    }

    #[test]
    #[should_panic]
    fn test_plane_intersect(){
        let plane = Plane::new(&CartessianND::new(vec![3f64; 2]), 0f64);

        let a = CartessianND::new(vec![-2f64; 2]);
        let movement = CartessianND::new(vec![0f64; 2]);
        plane.find_intersect(&a, &movement);
    }

    #[test]
    fn test_simpleplanepair(){
        // Float
        let planepair = SimplePlanePair::new(0, [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(planepair.pos[0], 0f64);
        assert_abs_diff_eq!(planepair.pos[1], 1f64);

        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([-1f64, 0.0])));
        assert_eq!(planepair.normal_at_unsafe(&a), Cartessian2D::new([-1f64, 0.0]));

        a[0] = 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([1f64, 0.0])));
        assert_eq!(planepair.normal_at_unsafe(&a), Cartessian2D::new([1f64, 0.0]));

        a[0] = -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0.5f64;
        let mut movement = Cartessian2D::new([-2f64; 2]);
        assert_eq!(planepair.find_intersect(&a, &movement), Some(Cartessian2D::new([0f64, 1.5f64])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([0f64, 1.5f64])));
        movement[0] = 2f64;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(Cartessian2D::new([1f64, 1.5f64])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([1f64, 1.5f64])));

        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![-1f64, 0.0])));
        assert_eq!(planepair.normal_at_unsafe(&a), CartessianND::new(vec![-1f64, 0.0]));

        a[0] = 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![1f64, 0.0])));
        assert_eq!(planepair.normal_at_unsafe(&a), CartessianND::new(vec![1f64, 0.0]));

        a[0] = -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0.5f64;
        let mut movement = CartessianND::new(vec![-2f64; 2]);
        assert_eq!(planepair.find_intersect(&a, &movement), Some(CartessianND::new(vec![0f64, 1.5f64])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(CartessianND::new(vec![0f64, 1.5f64])));
        movement[0] = 2f64;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(CartessianND::new(vec![1f64, 1.5f64])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(CartessianND::new(vec![1f64, 1.5f64])));

        // Integer
        let planepair = SimplePlanePair::new(0, [2, 0]).unwrap();
        assert_abs_diff_eq!(planepair.pos[0], 0);
        assert_abs_diff_eq!(planepair.pos[1], 2);

        let mut a = Cartessian2D::new([3; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 2;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([-1, 0])));
        assert_eq!(planepair.normal_at_unsafe(&a), Cartessian2D::new([-1, 0]));

        a[0] = 1;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([1, 0])));
        assert_eq!(planepair.normal_at_unsafe(&a), Cartessian2D::new([1, 0]));

        a[0] = -1;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0;
        let mut movement = Cartessian2D::new([-1, 0]);
        assert_eq!(planepair.find_intersect(&a, &movement), Some(a.clone()));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(a.clone()));

        a[0] = 2;
        movement[0] = 1;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(a.clone()));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(a.clone()));

        let mut a = CartessianND::new(vec![3; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 2;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![-1, 0])));
        assert_eq!(planepair.normal_at_unsafe(&a), CartessianND::new(vec![-1, 0]));

        a[0] = 1;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![1, 0])));
        assert_eq!(planepair.normal_at_unsafe(&a), CartessianND::new(vec![1, 0]));

        a[0] = -1;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0;
        let mut movement = CartessianND::new(vec![-1, 0]);
        assert_eq!(planepair.find_intersect(&a, &movement), Some(a.clone()));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(a.clone()));

        a[0] = 2;
        movement[0] = 1;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(a.clone()));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(a.clone()));


        assert_eq!(format!("{}", planepair), "SimplePlanePair between x[0] = 0 and 2");
        assert_eq!(format!("{:.3e}", planepair), "SimplePlanePair between x[0] = 0.000e0 and 2.000e0");
        assert_eq!(format!("{}", planepair.brief()), "0,0,2");
        assert_eq!(format!("{:.3e}", planepair.brief()), "0,0.000e0,2.000e0");
    }

    #[test]
    #[should_panic]
    fn test_planepair_panic(){
        let planepair = SimplePlanePair::new(0, [1f64, 0f64]).unwrap();
        let a = Cartessian2D::new([3f64; 2]);
        let movement = Cartessian2D::new([-2f64; 2]);
        planepair.find_intersect(&a, &movement);
    }

    #[test]
    fn test_planepair(){
        let planepair = PlanePair::new(Cartessian2D::new([1f64, 1f64]), [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(planepair.normal_vec, Cartessian2D::new([1.0 / 2f64.sqrt(); 2]));

        let mut a = Cartessian2D::new([1f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([-1.0 / 2f64.sqrt(); 2])));
        assert_eq!(planepair.normal_at_unsafe(&a), Cartessian2D::new([-1.0 / 2f64.sqrt(); 2]));

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a *= -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([1.0 / 2f64.sqrt(); 2])));
        assert_eq!(planepair.normal_at_unsafe(&a), Cartessian2D::new([1.0 / 2f64.sqrt(); 2]));

        let a = Cartessian2D::new([0f64, 0.5f64]);
        let mut movement = Cartessian2D::new([0f64, 1f64]);
        assert_eq!(planepair.find_intersect(&a, &movement), Some(Cartessian2D::new([0f64, 1f64])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([0f64, 1f64])));
        movement[1] = -1f64;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(Cartessian2D::new([0f64; 2])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([0f64; 2])));

        let planepair = PlanePair::new(CartessianND::new(vec![1f64, 1f64]), [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(planepair.normal_vec, CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]));

        let mut a = CartessianND::new(vec![1f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![-1.0 / 2f64.sqrt(); 2])));
        assert_eq!(planepair.normal_at_unsafe(&a), CartessianND::new(vec![-1.0 / 2f64.sqrt(); 2]));

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a *= -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])));
        assert_eq!(planepair.normal_at_unsafe(&a), CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]));

        let a = CartessianND::new(vec![0f64, 0.5f64]);
        let mut movement = CartessianND::new(vec![0f64, 1f64]);
        assert_eq!(planepair.find_intersect(&a, &movement), Some(CartessianND::new(vec![0f64, 1f64])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(CartessianND::new(vec![0f64, 1f64])));
        movement[1] = -1f64;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(CartessianND::new(vec![0f64; 2])));
        assert_eq!(planepair.find_intersect_unsafe(&a, &movement), Some(CartessianND::new(vec![0f64; 2])));

        assert_eq!(format!("{}", planepair), "PlanePair normal to (0.7071067811865475:0.7071067811865475) between constant 0 and 0.7071067811865475");
        assert_eq!(format!("{:.3e}", planepair), "PlanePair normal to (7.071e-1:7.071e-1) between constant 0.000e0 and 7.071e-1");
        assert_eq!(format!("{}", planepair.brief()), "0.7071067811865475:0.7071067811865475,0,0.7071067811865475");
        assert_eq!(format!("{:.3e}", planepair.brief()), "7.071e-1:7.071e-1,0.000e0,7.071e-1");
    }

    #[test]
    fn test_simplebox(){
        let simplebox = SimpleBox::cube_with_center(&Cartessian3D::new([0; 3]), 3).unwrap();
        for plane in &simplebox.planes{
            match plane.idx{
                0 => {
                    assert_eq!(plane, &SimplePlanePair::new(0, [-3, 3]).unwrap());
                },
                1 => {
                    assert_eq!(plane, &SimplePlanePair::new(1, [-3, 3]).unwrap());
                },
                2 => {
                    assert_eq!(plane, &SimplePlanePair::new(2, [-3, 3]).unwrap());
                },
                _ => {
                    unreachable!();
                }
            }
        }

        let mut a = Cartessian3D::new([0, 0, 0]);
        let mut vec = Cartessian3D::new([0, 0, 0]);

        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), None);

        a[0] = -3;
        vec[0] = 1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a[0] = 3;
        vec[0] = -1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a[0] = 5;
        assert_eq!(simplebox.check_inclusion(&a), false);
        assert_eq!(simplebox.normal_at(&a), None);

        a[0] = 0;
        a[1] = 3;
        vec[0] = 0;
        vec[1] = -1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a[1] = -3;
        vec[1] = 1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a[1] = 0; vec[1] = 0;
        a[2] = 3; vec[2] = -1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a[2] = -3; vec[2] = 1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a.coord = [3, 0, 0];
        vec.coord = [1, 0, 0];
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = -3; vec[0] = -1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = 0; vec[0] = 0;
        a[1] = 3; vec[1] = 1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = -3; vec[1] = -1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = 0; vec[1] = 0;
        a[2] = 3; vec[2] = 1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[2] = -3; vec[2] = -1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        assert_eq!(format!("{}", simplebox), "SimpleBox enveloped in plane pairs (0,-3,3), (1,-3,3), (2,-3,3)");
        assert_eq!(format!("{:.3e}", simplebox), "SimpleBox enveloped in plane pairs (0,-3.000e0,3.000e0), (1,-3.000e0,3.000e0), (2,-3.000e0,3.000e0)");
        assert_eq!(format!("{}", simplebox.brief()), "(0,-3,3),(1,-3,3),(2,-3,3)");
        assert_eq!(format!("{:.3e}", simplebox.brief()), "(0,-3.000e0,3.000e0),(1,-3.000e0,3.000e0),(2,-3.000e0,3.000e0)");
    }

    #[test]
    fn test_cube(){
        let cube = Cube::new(&Cartessian3D::new([0; 3]), 3);
        assert_eq!(&cube.center, &Cartessian3D::new([0; 3]));
        assert_eq!(&cube.radius, &3);

        let mut a = Cartessian3D::new([0, 0, 0]);
        let mut vec = Cartessian3D::new([0, 0, 0]);

        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), None);

        a[0] = -3;
        vec[0] = 1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[0] = 3;
        vec[0] = -1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[0] = 5;
        assert_eq!(cube.check_inclusion(&a), false);
        assert_eq!(cube.normal_at(&a), None);

        a[0] = 0;
        a[1] = 3;
        vec[0] = 0;
        vec[1] = -1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[1] = -3;
        vec[1] = 1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[1] = 0; vec[1] = 0;
        a[2] = 3; vec[2] = -1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[2] = -3; vec[2] = 1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a.coord = [3, 0, 0];
        vec.coord = [1, 0, 0];
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = -3; vec[0] = -1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = 0; vec[0] = 0;
        a[1] = 3; vec[1] = 1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = -3; vec[1] = -1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = 0; vec[1] = 0;
        a[2] = 3; vec[2] = 1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[2] = -3; vec[2] = -1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        assert_eq!(format!("{}", cube), "Cube has center at (0:0:0) with radius 3");
        assert_eq!(format!("{:.3e}", cube), "Cube has center at (0.000e0:0.000e0:0.000e0) with radius 3.000e0");
        assert_eq!(format!("{}", cube.brief()), "0:0:0,3");
        assert_eq!(format!("{:.3e}", cube.brief()), "0.000e0:0.000e0:0.000e0,3.000e0");
    }

    #[test]
    fn test_parallelogram(){
        let vec1 = Cartessian2D::new([0f64, 1f64]);
        let mut vec2 = Cartessian2D::new([1f64, -1f64]);

        let plane1 = PlanePair::new(vec1.clone(), [0f64, 1f64]).unwrap();
        let plane2 = PlanePair::new(vec2.clone(), [0f64, 1f64]).unwrap();
        let parallelogram = Parallelogram::new([plane1.clone(), plane2.clone()]).unwrap();

        let mut a = Cartessian2D::new([1f64, 0.5f64]);
        vec2 /= 2f64.sqrt();

        assert_eq!(parallelogram.check_inclusion(&a), true);
        assert_eq!(parallelogram.normal_at(&a), None);

        a[0] = 0.5f64;
        assert_eq!(parallelogram.check_inclusion(&a), true);
        assert_eq!(parallelogram.normal_at(&a), Some(vec2.clone()));
        assert_eq!(parallelogram.normal_at_unsafe(&a), vec2.clone());

        a[0] = 1.5f64;
        assert_eq!(parallelogram.check_inclusion(&a), true);
        assert_eq!(parallelogram.normal_at(&a), Some(-vec2.clone()));
        assert_eq!(parallelogram.normal_at_unsafe(&a), -vec2.clone());

        a[0] = 5f64;
        assert_eq!(parallelogram.check_inclusion(&a), false);
        assert_eq!(parallelogram.normal_at(&a), None);

        a[0] = 1f64;
        a[1] = 1f64;
        assert_eq!(parallelogram.check_inclusion(&a), true);
        assert_eq!(parallelogram.normal_at(&a), Some(-vec1.clone()));
        assert_eq!(parallelogram.normal_at_unsafe(&a), -vec1.clone());

        a[1] = 0f64;
        assert_eq!(parallelogram.check_inclusion(&a), true);
        assert_eq!(parallelogram.normal_at(&a), Some(vec1.clone()));
        assert_eq!(parallelogram.normal_at_unsafe(&a), vec1.clone());

        a[1] = -1f64;
        assert_eq!(parallelogram.check_inclusion(&a), false);
        assert_eq!(parallelogram.normal_at(&a), None);

        a.coord = [1f64, 0.5f64];
        let mut movement = Cartessian2D::new([1f64, 1f64]);
        assert_eq!(parallelogram.find_intersect(&a, &movement), Some(Cartessian2D::new([1.5f64, 1f64])));
        assert_eq!(parallelogram.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([1.5f64, 1f64])));

        movement.coord = [-1f64; 2];
        assert_eq!(parallelogram.find_intersect(&a, &movement), Some(Cartessian2D::new([0.5f64, 0f64])));
        assert_eq!(parallelogram.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([0.5f64, 0f64])));

        movement.coord = [2f64, 0f64];
        assert_eq!(parallelogram.find_intersect(&a, &movement), Some(Cartessian2D::new([1.5f64, 0.5f64])));
        assert_eq!(parallelogram.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([1.5f64, 0.5f64])));

        movement.coord = [-2f64, 0f64];
        assert_eq!(parallelogram.find_intersect(&a, &movement), Some(Cartessian2D::new([0.5f64, 0.5f64])));
        assert_eq!(parallelogram.find_intersect_unsafe(&a, &movement), Some(Cartessian2D::new([0.5f64, 0.5f64])));

        assert_eq!(format!("{}", parallelogram), "Parallelogram enveloped in plane pairs (0:1,0,1), (0.7071067811865475:-0.7071067811865475,0,0.7071067811865475)");
        assert_eq!(format!("{:.3e}", parallelogram), "Parallelogram enveloped in plane pairs (0.000e0:1.000e0,0.000e0,1.000e0), (7.071e-1:-7.071e-1,0.000e0,7.071e-1)");
        assert_eq!(format!("{}", parallelogram.brief()), "(0:1,0,1),(0.7071067811865475:-0.7071067811865475,0,0.7071067811865475)");
        assert_eq!(format!("{:.3e}", parallelogram.brief()), "(0.000e0:1.000e0,0.000e0,1.000e0),(7.071e-1:-7.071e-1,0.000e0,7.071e-1)");
    }
}
