// use crate::boundary::BoundaryCondition;
// use crate::boundary::AfterMove;
use crate::argument::CommandBuilder;
use crate::boundary::NonPeriodic;
use crate::boundary::Periodic;
use crate::boundary::{FloatBoundary, IntBoundary};
use crate::vector::basic::Map;
use crate::vector::product::Norm;
use crate::{prelude::Error, vector::Dim};
use approx::AbsDiffEq;
use clap::Arg;
use clap::ArgMatches;
use std::convert::TryInto;
use std::fmt::Debug;
use std::fmt::Formatter;
use std::fmt::{self, Display};
use std::ops::Div;
use std::ops::Rem;

use crate::vector::product::Dot;
use crate::vector::Scalar;
use crate::vector::Vector;
use crate::vector::{Cartessian, CartessianND};
use serde::{Deserialize, Serialize};
use serde_json::from_str;
use std::ops::Neg;

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Direction {
    Positive,
    Negative,
}

impl Display for Direction {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            &Direction::Positive => {
                write!(f, "Positive")
            }
            &Direction::Negative => {
                write!(f, "Negative")
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SimplePlane<T> {
    idx: usize,
    pos: T,
    dir_in: Direction,
}

impl<T> SimplePlane<T> {
    pub fn new(idx: usize, pos: T, dir_in: Direction) -> Self {
        Self { idx, pos, dir_in }
    }
}

impl<T> NonPeriodic for SimplePlane<T> {}

macro_rules! impl_arg_simpleplane {
    ($ty : ident) => {
        impl<'h> CommandBuilder<'h, 3> for SimplePlane<$ty> {
            const SUBCOMMAND: &'h str = "SimplePlane";

            fn args() -> [Arg<'h>; 3] {
                [
                    Arg::new("idx")
                        .short('i')
                        .long("idx")
                        .value_name("IDX")
                        .takes_value(true)
                        .value_parser(clap::value_parser!(usize))
                        .help("Index of plane position"),
                    Arg::new("pos")
                        .short('p')
                        .long("pos")
                        .value_name("POS")
                        .takes_value(true)
                        .value_parser(clap::value_parser!($ty))
                        .help("Position of plane"),
                    Arg::new("direction")
                        .short('d')
                        .long("dir")
                        .value_name("DIR")
                        .takes_value(true)
                        .value_parser(clap::value_parser!(usize))
                        .help("Direction to interior. 0 : Negative, 1 : Positive"),
                ]
            }
        }

        impl From<&ArgMatches> for SimplePlane<$ty> {
            fn from(m: &ArgMatches) -> Self {
                let idx: usize = *m.get_one::<usize>("idx").unwrap();
                let pos: $ty = *m.get_one::<$ty>("pos").unwrap();
                let dir: Direction = match m.get_one::<usize>("direction") {
                    Some(&0) => Direction::Negative,
                    Some(&1) => Direction::Positive,
                    Some(&_) | None => {
                        panic!("Invalid argument input for direction. Only 0 and 1 are valid.");
                    }
                };
                SimplePlane::<$ty>::new(idx, pos, dir)
            }
        }
    };
}

impl_arg_simpleplane!(f32);
impl_arg_simpleplane!(f64);
impl_arg_simpleplane!(i8);
impl_arg_simpleplane!(i16);
impl_arg_simpleplane!(i32);
impl_arg_simpleplane!(i64);
impl_arg_simpleplane!(i128);
impl_arg_simpleplane!(isize);

macro_rules! impl_float_simpleplane {
    ($ty : ident) => {
        impl<const N : usize> FloatBoundary<Cartessian<$ty, N>> for SimplePlane<$ty>{
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

            fn ratio_to_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<$ty>{
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
                        return Some(t);
                    },
                    Direction::Negative => {
                        let t = (self.pos - pos[self.idx]) / movement[self.idx];
                        return Some(t);
                    },
                }
            }

            fn ratio_to_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<$ty>{
                let destination = pos + movement;
                if self.check_inclusion(&destination){
                    return None;
                }

                match self.dir_in{
                    Direction::Positive => {
                        let t = (self.pos - pos[self.idx]) / movement[self.idx];
                        return Some(t);
                    },
                    Direction::Negative => {
                        let t = (self.pos - pos[self.idx]) / movement[self.idx];
                        return Some(t);
                    },
                }
            }
        }

        impl FloatBoundary<CartessianND<$ty>> for SimplePlane<$ty>{
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

            fn ratio_to_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<$ty>{
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
                            return Some(t);
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            let t = (self.pos - pos[self.idx]) / movement[self.idx];
                            return Some(t);
                        }
                    },
                }
            }

            fn ratio_to_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<$ty>{
                match self.dir_in{
                    Direction::Positive => {
                        if pos[self.idx] + movement[self.idx] >= self.pos{
                            return None;
                        } else {
                            let t = (self.pos - pos[self.idx]) / movement[self.idx];
                            return Some(t);
                        }
                    },
                    Direction::Negative => {
                        if pos[self.idx] + movement[self.idx] <= self.pos{
                            return None;
                        } else {
                            let t = (self.pos - pos[self.idx]) / movement[self.idx];
                            return Some(t);
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
        impl<const N : usize> IntBoundary<Cartessian<$ty, N>> for SimplePlane<$ty>{
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

        impl IntBoundary<CartessianND<$ty>> for SimplePlane<$ty>{
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

// impl<T: Display> Display for SimplePlane<T> {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         match self.dir_in {
//             Direction::Positive => write!(
//                 f,
//                 "SimplePlane at x[{}] = {} with positive to inside",
//                 self.idx, self.pos
//             ),
//             Direction::Negative => write!(
//                 f,
//                 "SimplePlane at x[{}] = {} with negative to inside",
//                 self.idx, self.pos
//             ),
//         }
//     }
// }

// impl<T: LowerExp> LowerExp for SimplePlane<T> {
//     fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
//         match self.dir_in {
//             Direction::Positive => {
//                 write!(f, "SimplePlane at x[{}] = ", self.idx)?;
//                 LowerExp::fmt(&self.pos, f)?;
//                 write!(f, " with positive to inside")
//             }
//             Direction::Negative => {
//                 write!(f, "SimplePlane at x[{}] = ", self.idx)?;
//                 LowerExp::fmt(&self.pos, f)?;
//                 write!(f, " with negative to inside")
//             }
//         }
//     }
// }

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Plane<V: Vector> {
    normal_vec: V, // direction of normal vec : inside
    constant: <V as Vector>::Item,
}

impl<V> Plane<V>
where
    V: Vector
        + Clone
        + Norm
        + Div<<V as Norm>::Output, Output = V>
        + Dot<V, Output = <V as Vector>::Item>,
    <V as Vector>::Item: Clone + Div<<V as Norm>::Output, Output = <V as Vector>::Item> + Debug,
    <V as Norm>::Output: Copy,
{
    pub fn new(normal_vec: V, constant: <V as Vector>::Item) -> Self {
        let n = normal_vec.norm_l2();
        Self {
            normal_vec: normal_vec.clone() / n,
            constant: constant / n,
        }
    }

    pub fn from_point(normal_vec: V, point: V) -> Self {
        let n = normal_vec.norm_l2();
        let vec = normal_vec.clone() / n;
        let c = vec.dot(point);
        Self {
            normal_vec: vec,
            constant: c / n,
        }
    }
}

impl<V: Vector> NonPeriodic for Plane<V> {}

macro_rules! impl_clap_plane {
    (fixed, $ty : ident) => {
        impl<'h, const N: usize> CommandBuilder<'h, 2> for Plane<Cartessian<$ty, N>> {
            const SUBCOMMAND: &'h str = "Plane";

            fn args() -> [Arg<'h>; 2] {
                [
                    Arg::new("normal_vec")
                        .short('n')
                        .long("normal")
                        .value_name("NORMAL")
                        .takes_value(true)
                        .help("Normal vector of plane."),
                    Arg::new("constant")
                        .short('c')
                        .long("const")
                        .value_name("CONST")
                        .takes_value(true)
                        .value_parser(clap::value_parser!($ty))
                        .help("Constant of plane."),
                ]
            }
        }

        impl<const N: usize> From<&ArgMatches> for Plane<Cartessian<$ty, N>> {
            fn from(m: &ArgMatches) -> Self {
                let normal: Cartessian<$ty, N> =
                    m.get_one::<String>("normal_vec").unwrap().parse().unwrap();
                let constant: $ty = *m.get_one::<$ty>("constant").unwrap();
                Plane::<Cartessian<$ty, N>>::new(normal, constant)
            }
        }
    };
    (nd, $ty : ident) => {
        impl<'h> CommandBuilder<'h, 2> for Plane<CartessianND<$ty>> {
            const SUBCOMMAND: &'h str = "Plane";

            fn args() -> [Arg<'h>; 2] {
                [
                    Arg::new("normal_vec")
                        .short('n')
                        .long("normal")
                        .value_name("NORMAL")
                        .takes_value(true)
                        .help("Normal vector of plane."),
                    Arg::new("constant")
                        .short('c')
                        .long("const")
                        .value_name("CONST")
                        .takes_value(true)
                        .value_parser(clap::value_parser!($ty))
                        .help("Constant of plane."),
                ]
            }
        }

        impl From<&ArgMatches> for Plane<CartessianND<$ty>> {
            fn from(m: &ArgMatches) -> Self {
                let normal: CartessianND<$ty> =
                    m.get_one::<String>("normal_vec").unwrap().parse().unwrap();
                let constant: $ty = *m.get_one::<$ty>("constant").unwrap();
                Plane::<CartessianND<$ty>>::new(normal, constant)
            }
        }
    };
    ($ty : ident) => {
        impl_clap_plane!(fixed, $ty);
        impl_clap_plane!(nd, $ty);
    };
}

impl_clap_plane!(f32);
impl_clap_plane!(f64);

macro_rules! impl_float_plane {
    ($ty : ident) => {
        impl<const N: usize> FloatBoundary<Cartessian<$ty, N>> for Plane<Cartessian<$ty, N>> {
            fn check_inclusion(&self, pos: &Cartessian<$ty, N>) -> bool {
                let d = self.normal_vec.dot(pos);

                d >= self.constant
            }

            fn normal_at(&self, pos: &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let d = self.normal_vec.dot(pos);
                if self
                    .constant
                    .abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon())
                {
                    return Some(self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, _pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                self.normal_vec.clone()
            }

            fn find_intersect(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let destination = pos + movement;
                if self.check_inclusion(&destination) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(pos + movement * t);
            }

            fn find_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                let destination = pos + movement;
                if self.check_inclusion(&destination) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(pos + movement * t);
            }

            fn ratio_to_intersect(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<$ty> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let destination = pos + movement;
                if self.check_inclusion(&destination) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(t);
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<$ty> {
                let destination = pos + movement;
                if self.check_inclusion(&destination) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(t);
            }
        }

        impl FloatBoundary<CartessianND<$ty>> for Plane<CartessianND<$ty>> {
            fn check_inclusion(&self, pos: &CartessianND<$ty>) -> bool {
                let d = self.normal_vec.dot(pos);

                d >= self.constant
            }

            fn normal_at(&self, pos: &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                let d = self.normal_vec.dot(pos);
                if self
                    .constant
                    .abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon())
                {
                    return Some(self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, _pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                self.normal_vec.clone()
            }

            fn find_intersect(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let mut result = pos.clone();
                result.zip_mut_with(movement, move |x, y| *x = *x + *y);
                if self.check_inclusion(&result) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                result.clone_from(pos);
                result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                return Some(result);
            }

            fn find_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                let mut result = pos.clone();
                result.zip_mut_with(movement, move |x, y| *x = *x + *y);
                if self.check_inclusion(&result) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                result.clone_from(pos);
                result.zip_mut_with(movement, move |x, y| *x = *x + *y * t);
                return Some(result);
            }

            fn ratio_to_intersect(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<$ty> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let mut result = pos.clone();
                result.zip_mut_with(movement, move |x, y| *x = *x + *y);
                if self.check_inclusion(&result) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(t);
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<$ty> {
                let mut result = pos.clone();
                result.zip_mut_with(movement, move |x, y| *x = *x + *y);
                if self.check_inclusion(&result) {
                    return None;
                }

                let t = -(self.normal_vec.dot(pos) + self.constant) / self.normal_vec.dot(movement);
                return Some(t);
            }
        }
    };
}

impl_float_plane!(f32);
impl_float_plane!(f64);

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SimplePlanePair<T> {
    idx: usize,
    pos: [T; 2],
}

impl<T: PartialOrd + Copy + AbsDiffEq> SimplePlanePair<T> {
    pub fn new(idx: usize, pos: [T; 2]) -> Option<Self> {
        if pos[0].abs_diff_eq(&pos[1], <T as AbsDiffEq>::default_epsilon()) {
            return None;
        }

        let (a, b) = match pos[0] > pos[1] {
            true => (pos[1], pos[0]),
            false => (pos[0], pos[1]),
        };

        Some(Self { idx, pos: [a, b] })
    }
}

macro_rules! impl_arg_simpleplanepair {
    ($ty : ident) => {
        impl<'h> CommandBuilder<'h, 2> for SimplePlanePair<$ty> {
            const SUBCOMMAND: &'h str = "SimplePlanePair";

            fn args() -> [Arg<'h>; 2] {
                [
                    Arg::new("idx")
                        .short('i')
                        .long("idx")
                        .value_name("IDX")
                        .takes_value(true)
                        .value_parser(clap::value_parser!(usize))
                        .help("Index of plane position"),
                    Arg::new("pos")
                        .short('p')
                        .long("pos")
                        .value_name("POS")
                        .takes_value(true)
                        .help("Positions of plane. format : [c1,c2]"),
                ]
            }
        }

        impl From<&ArgMatches> for SimplePlanePair<$ty> {
            fn from(m: &ArgMatches) -> Self {
                let idx: usize = *m.get_one::<usize>("idx").unwrap();
                let pos: [$ty; 2] = from_str(m.get_one::<String>("pos").unwrap()).unwrap();
                SimplePlanePair::<$ty>::new(idx, pos).unwrap()
            }
        }
    };
}

impl_arg_simpleplanepair!(f32);
impl_arg_simpleplanepair!(f64);
impl_arg_simpleplanepair!(i8);
impl_arg_simpleplanepair!(i16);
impl_arg_simpleplanepair!(i32);
impl_arg_simpleplanepair!(i64);
impl_arg_simpleplanepair!(isize);

macro_rules! impl_float_simpleplanepair {
    ($ty : ident) => {
        impl<const N : usize> FloatBoundary<Cartessian<$ty, N>> for SimplePlanePair<$ty>{
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

            fn ratio_to_intersect(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<$ty> {
                if self.idx >= N{
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else {
                    return None;
                }
            }

            fn ratio_to_intersect_unsafe(&self, pos : &Cartessian<$ty, N>, movement : &Cartessian<$ty, N>) -> Option<$ty> {
                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else {
                    return None;
                }
            }
        }

        impl<const N : usize> Periodic<Cartessian<$ty, N>> for SimplePlanePair<$ty>{
            fn find_pair(&self, pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                if self.idx >= N {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut result = pos.clone();
                if result[self.idx] < self.pos[0]{
                    result[self.idx] += self.pos[1] - self.pos[0];
                } else if self.pos[1] < result[self.idx]{
                    result[self.idx] -= self.pos[1] - self.pos[0];
                }

                return result;
            }

            fn find_pair_mut(&self, pos : &mut Cartessian<$ty, N>) {
                if self.idx >= N {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if pos[self.idx] < self.pos[0]{
                    pos[self.idx] += self.pos[1] - self.pos[0];
                } else if self.pos[1] < pos[self.idx]{
                    pos[self.idx] -= self.pos[1] - self.pos[0];
                }
            }
        }

        impl FloatBoundary<CartessianND<$ty>> for SimplePlanePair<$ty>{
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

            fn ratio_to_intersect(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<$ty> {
                if self.idx >= pos.dim(){
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                } else if !self.check_inclusion(pos){
                    panic!("State cannot live outside of system boundary. Move from outside occurs");
                }

                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else {
                    return None;
                }
            }

            fn ratio_to_intersect_unsafe(&self, pos : &CartessianND<$ty>, movement : &CartessianND<$ty>) -> Option<$ty> {
                let destination = pos[self.idx] + movement[self.idx];
                if destination < self.pos[0]{
                    let t = (self.pos[0] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else if self.pos[1] < destination{
                    let t = (self.pos[1] - pos[self.idx]) / movement[self.idx];
                    return Some(t);
                } else {
                    return None;
                }
            }
        }

        impl Periodic<CartessianND<$ty>> for SimplePlanePair<$ty>{
            fn find_pair(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty> {
                if self.idx >= pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut result = pos.clone();
                if result[self.idx] < self.pos[0]{
                    result[self.idx] += self.pos[1] - self.pos[0];
                } else if self.pos[1] < result[self.idx]{
                    result[self.idx] -= self.pos[1] - self.pos[0];
                }

                return result;
            }

            fn find_pair_mut(&self, pos : &mut CartessianND<$ty>) {
                if self.idx >= pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if pos[self.idx] < self.pos[0]{
                    pos[self.idx] += self.pos[1] - self.pos[0];
                } else if self.pos[1] < pos[self.idx]{
                    pos[self.idx] -= self.pos[1] - self.pos[0];
                }
            }
        }
    };
}

impl_float_simpleplanepair!(f32);
impl_float_simpleplanepair!(f64);

macro_rules! impl_int_simpleplanepair {
    ($ty : ident) => {
        impl<const N : usize> IntBoundary<Cartessian<$ty, N>> for SimplePlanePair<$ty>{
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

        impl<const N : usize> Periodic<Cartessian<$ty, N>> for SimplePlanePair<$ty>{
            fn find_pair(&self, pos : &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                if self.idx >= N {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut result = pos.clone();
                if result[self.idx] < self.pos[0]{
                    result[self.idx] += self.pos[1] - self.pos[0] + 1;
                } else if self.pos[1] < result[self.idx]{
                    result[self.idx] -= self.pos[1] - self.pos[0] + 1;
                }

                return result;
            }

            fn find_pair_mut(&self, pos : &mut Cartessian<$ty, N>) {
                if self.idx >= N {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if pos[self.idx] < self.pos[0]{
                    pos[self.idx] += self.pos[1] - self.pos[0] + 1;
                } else if self.pos[1] < pos[self.idx]{
                    pos[self.idx] -= self.pos[1] - self.pos[0] + 1;
                }
            }
        }

        impl IntBoundary<CartessianND<$ty>> for SimplePlanePair<$ty>{
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

        impl Periodic<CartessianND<$ty>> for SimplePlanePair<$ty>{
            fn find_pair(&self, pos : &CartessianND<$ty>) -> CartessianND<$ty> {
                if self.idx >= pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                let mut result = pos.clone();
                if result[self.idx] < self.pos[0]{
                    result[self.idx] += self.pos[1] - self.pos[0] + 1;
                } else if self.pos[1] < result[self.idx]{
                    result[self.idx] -= self.pos[1] - self.pos[0] + 1;
                }

                return result;
            }

            fn find_pair_mut(&self, pos : &mut CartessianND<$ty>) {
                if self.idx >= pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
                }

                if pos[self.idx] < self.pos[0]{
                    pos[self.idx] += self.pos[1] - self.pos[0] + 1;
                } else if self.pos[1] < pos[self.idx]{
                    pos[self.idx] -= self.pos[1] - self.pos[0] + 1;
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

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PlanePair<V: Vector> {
    normal_vec: V,
    constant: [<V as Vector>::Item; 2],
}

impl<V> PlanePair<V>
where
    V: Vector + Clone + Norm + Div<<V as Norm>::Output, Output = V>,
    <V as Vector>::Item:
        PartialOrd + AbsDiffEq + Div<<V as Norm>::Output, Output = <V as Vector>::Item>,
    <V as Norm>::Output: Copy,
{
    pub fn new(normal_vec: V, constant: [<V as Vector>::Item; 2]) -> Option<Self> {
        if constant[0].abs_diff_eq(
            &constant[1],
            <<V as Vector>::Item as AbsDiffEq>::default_epsilon(),
        ) {
            return None;
        }

        let n = normal_vec.norm_l2();

        let (a, b) = match constant[0] > constant[1] {
            true => (constant[1], constant[0]),
            false => (constant[0], constant[1]),
        };
        Some(Self {
            normal_vec: normal_vec.clone() / n,
            constant: [a / n, b / n],
        })
    }
}

macro_rules! impl_planepair_clap {
    (fixed, $ty : ident) => {
        impl<'h, const N: usize> CommandBuilder<'h, 2> for PlanePair<Cartessian<$ty, N>> {
            const SUBCOMMAND: &'h str = "PlanePair";

            fn args() -> [Arg<'h>; 2] {
                [
                    Arg::new("normal_vec")
                        .short('n')
                        .long("normal")
                        .value_name("NORMAL")
                        .takes_value(true)
                        .help("Normal vector of plane."),
                    Arg::new("constant")
                        .short('c')
                        .long("const")
                        .value_name("CONST")
                        .takes_value(true)
                        .help("Constant of plane. format : [c1,c2]"),
                ]
            }
        }

        impl<const N: usize> From<&ArgMatches> for PlanePair<Cartessian<$ty, N>> {
            fn from(m: &ArgMatches) -> Self {
                let normal: Cartessian<$ty, N> =
                    m.get_one::<String>("normal_vec").unwrap().parse().unwrap();
                let constant: [$ty; 2] = from_str(m.get_one::<String>("constant").unwrap()).unwrap();
                PlanePair::<Cartessian<$ty, N>>::new(normal, constant).unwrap()
            }
        }
    };
    (nd, $ty : ident) => {
        impl<'h> CommandBuilder<'h, 2> for PlanePair<CartessianND<$ty>> {
            const SUBCOMMAND: &'h str = "PlanePair";

            fn args() -> [Arg<'h>; 2] {
                [
                    Arg::new("normal_vec")
                        .short('n')
                        .long("normal")
                        .value_name("NORMAL")
                        .takes_value(true)
                        .help("Normal vector of plane."),
                    Arg::new("constant")
                        .short('c')
                        .long("const")
                        .value_name("CONST")
                        .takes_value(true)
                        .help("Constant of plane. format : [c1,c2]"),
                ]
            }
        }

        impl From<&ArgMatches> for PlanePair<CartessianND<$ty>> {
            fn from(m: &ArgMatches) -> Self {
                let normal: CartessianND<$ty> =
                    m.get_one::<String>("normal_vec").unwrap().parse().unwrap();
                let constant: [$ty; 2] = from_str(m.get_one::<String>("constant").unwrap()).unwrap();
                PlanePair::<CartessianND<$ty>>::new(normal, constant).unwrap()
            }
        }
    };
    ($ty : ident) => {
        impl_planepair_clap!(fixed, $ty);
        impl_planepair_clap!(nd, $ty);
    };
}

impl_planepair_clap!(f32);
impl_planepair_clap!(f64);

macro_rules! impl_float_planepair {
    ($ty : ident) => {
        impl<const N: usize> FloatBoundary<Cartessian<$ty, N>> for PlanePair<Cartessian<$ty, N>> {
            fn check_inclusion(&self, pos: &Cartessian<$ty, N>) -> bool {
                let d = self.normal_vec.dot(pos);

                self.constant[0] <= d && d <= self.constant[1]
            }

            fn normal_at(&self, pos: &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()) {
                    return Some(self.normal_vec.clone());
                } else if self.constant[1].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()) {
                    return Some(-self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()) {
                    return self.normal_vec.clone();
                } else {
                    return -self.normal_vec.clone();
                }
            }

            fn find_intersect(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let d = self.normal_vec.dot(pos + movement);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                let d = self.normal_vec.dot(pos + movement);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(pos + movement * t);
                } else {
                    return None;
                }
            }

            fn ratio_to_intersect(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<$ty> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let d = self.normal_vec.dot(pos + movement);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else {
                    return None;
                }
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<$ty> {
                let d = self.normal_vec.dot(pos + movement);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else {
                    return None;
                }
            }
        }

        impl<const N: usize> Periodic<Cartessian<$ty, N>> for PlanePair<Cartessian<$ty, N>> {
            fn find_pair(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let d = self.normal_vec.dot(pos);
                let mut res = pos.clone();

                if d < self.constant[0] {
                    let dx = self.constant[1] - self.constant[0];
                    res += dx * &self.normal_vec;
                } else if self.constant[1] < d {
                    let dx = self.constant[1] - self.constant[0];
                    res -= dx * &self.normal_vec;
                }

                return res;
            }

            fn find_pair_mut(&self, pos: &mut Cartessian<$ty, N>) {
                let d = self.normal_vec.dot(&*pos);

                if d < self.constant[0] {
                    let dx = self.constant[1] - self.constant[0];
                    *pos += dx * &self.normal_vec;
                } else if self.constant[1] < d {
                    let dx = self.constant[1] - self.constant[0];
                    *pos -= dx * &self.normal_vec;
                }
            }
        }

        impl FloatBoundary<CartessianND<$ty>> for PlanePair<CartessianND<$ty>> {
            fn check_inclusion(&self, pos: &CartessianND<$ty>) -> bool {
                let d = self.normal_vec.dot(pos);

                self.constant[0] <= d && d <= self.constant[1]
            }

            fn normal_at(&self, pos: &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()) {
                    return Some(self.normal_vec.clone());
                } else if self.constant[1].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()) {
                    return Some(-self.normal_vec.clone());
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let d = self.normal_vec.dot(pos);

                if self.constant[0].abs_diff_eq(&d, <$ty as AbsDiffEq>::default_epsilon()) {
                    return self.normal_vec.clone();
                } else {
                    return -self.normal_vec.clone();
                }
            }

            fn find_intersect(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let mut result = pos.clone();
                result.zip_mut_with(movement, |x, y| *x = *x + *y);
                let d = self.normal_vec.dot(&result);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else {
                    return None;
                }
            }

            fn find_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                let mut result = pos.clone();
                result.zip_mut_with(movement, |x, y| *x = *x + *y);
                let d = self.normal_vec.dot(&result);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    result.clone_from(pos);
                    result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                    return Some(result);
                } else {
                    return None;
                }
            }

            fn ratio_to_intersect(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<$ty> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let mut result = pos.clone();
                result.zip_mut_with(movement, |x, y| *x = *x + *y);
                let d = self.normal_vec.dot(&result);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else {
                    return None;
                }
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<$ty> {
                let mut result = pos.clone();
                result.zip_mut_with(movement, |x, y| *x = *x + *y);
                let d = self.normal_vec.dot(&result);
                if d < self.constant[0] {
                    let t = (self.constant[0] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else if self.constant[1] < d {
                    let t = (self.constant[1] - self.normal_vec.dot(pos))
                        / self.normal_vec.dot(movement);
                    return Some(t);
                } else {
                    return None;
                }
            }
        }

        impl Periodic<CartessianND<$ty>> for PlanePair<CartessianND<$ty>> {
            fn find_pair(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let d = self.normal_vec.dot(pos);
                let mut res = pos.clone();

                if d < self.constant[0] {
                    let dx = self.constant[1] - self.constant[0];
                    res.zip_mut_with(&self.normal_vec, |x, y| *x = *x + *y * dx);
                } else if self.constant[1] < d {
                    let dx = self.constant[1] - self.constant[0];
                    res.zip_mut_with(&self.normal_vec, |x, y| *x = *x - *y * dx);
                }

                return res;
            }

            fn find_pair_mut(&self, pos: &mut CartessianND<$ty>) {
                let d = self.normal_vec.dot(&*pos);

                if d < self.constant[0] {
                    let dx = self.constant[1] - self.constant[0];
                    pos.zip_mut_with(&self.normal_vec, |x, y| *x = *x + *y * dx);
                } else if self.constant[1] < d {
                    let dx = self.constant[1] - self.constant[0];
                    pos.zip_mut_with(&self.normal_vec, |x, y| *x = *x - *y * dx);
                }
            }
        }
    };
}

impl_float_planepair!(f32);
impl_float_planepair!(f64);

#[derive(Clone, Debug, PartialEq)]
pub struct SimpleBox<T, const N: usize> {
    pub planes: [SimplePlanePair<T>; N],
}

macro_rules! impl_clap_simplebox {
    ($ty : ident, $n : expr) => {
        impl<'h> CommandBuilder<'h, 1> for SimpleBox<$ty, $n> {
            const SUBCOMMAND: &'h str = "SimpleBox";

            fn args() -> [Arg<'h>; 1] {
                [Arg::new("pairs")
                    .short('p')
                    .long("pairs")
                    .value_name("PAIRS")
                    .takes_value(true)
                    .help("Constant pairs. format : [[1.0, 2.0], [3.0, 4.0]]"),]
            }
        }

        impl From<&ArgMatches> for SimpleBox<$ty, $n> {
            fn from(m: &ArgMatches) -> Self {
                let pairs: [[$ty; 2]; $n] = from_str(m.get_one::<String>("pairs").unwrap()).unwrap();
                SimpleBox::<$ty, $n>::from_pairs(pairs)
            }
        }
    };
    ($ty : ident) => {
        impl_clap_simplebox!($ty, 1);
        impl_clap_simplebox!($ty, 2);
        impl_clap_simplebox!($ty, 3);
        impl_clap_simplebox!($ty, 4);
        impl_clap_simplebox!($ty, 5);
        impl_clap_simplebox!($ty, 6);
        impl_clap_simplebox!($ty, 7);
        impl_clap_simplebox!($ty, 8);
        impl_clap_simplebox!($ty, 9);
        impl_clap_simplebox!($ty, 10);
    };
}

impl_clap_simplebox!(f32);
impl_clap_simplebox!(f64);
impl_clap_simplebox!(i8);
impl_clap_simplebox!(i16);
impl_clap_simplebox!(i32);
impl_clap_simplebox!(i64);
impl_clap_simplebox!(isize);


impl<T, const N: usize> SimpleBox<T, N> {
    pub fn new(planes: [SimplePlanePair<T>; N]) -> Result<Self, Error>
    where
        T: Clone,
    {
        use itertools::Itertools;

        let count = planes.iter().map(|x| x.idx).unique().count();
        if count != N {
            return Err(Error::make_error_syntax(
                crate::prelude::ErrorCode::InvalidArgumentInput,
            ));
        }

        Ok(Self {
            planes: planes.clone(),
        })
    }

    pub fn from_pairs(consts: [[T; 2]; N]) -> Self
    where
        T: PartialOrd + Copy + AbsDiffEq + Debug,
    {
        let planes: [SimplePlanePair<T>; N] = consts
            .iter()
            .enumerate()
            .map(|(idx, x)| SimplePlanePair::new(idx, *x).unwrap())
            .collect::<Vec<SimplePlanePair<T>>>()
            .try_into()
            .unwrap();
        Self { planes }
    }

    pub fn cube_with_center<'a, V>(center: &'a V, half: T) -> Result<Self, Error>
    where
        V: Vector<Item = T> + Dim<N>,
        &'a V: IntoIterator<Item = &'a T>,
        T: Scalar + PartialOrd + AbsDiffEq + Rem<Output = T>,
    {
        if center.dim() != N || half < T::zero() {
            return Err(Error::make_error_syntax(
                crate::prelude::ErrorCode::InvalidDimension,
            ));
        }

        let planes: [SimplePlanePair<T>; N] = (0..N)
            .zip(center.into_iter())
            .map(|(idx, y)| SimplePlanePair::new(idx, [*y - half, *y + half]).unwrap())
            .collect::<Vec<SimplePlanePair<T>>>()
            .try_into()
            .unwrap();
        Ok(Self { planes })
    }
}

impl<T, const N: usize> FloatBoundary<Cartessian<T, N>> for SimpleBox<T, N>
where
    T: Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
    SimplePlanePair<T>: FloatBoundary<Cartessian<T, N>>,
{
    fn check_inclusion(&self, pos: &Cartessian<T, N>) -> bool {
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos: &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos: &Cartessian<T, N>) -> Cartessian<T, N> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return v,
                None => {}
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<Cartessian<T, N>> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn find_intersect_unsafe(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<Cartessian<T, N>> {
        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn ratio_to_intersect(&self, pos: &Cartessian<T, N>, movement: &Cartessian<T, N>) -> Option<T> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.ratio_to_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn ratio_to_intersect_unsafe(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<T> {
        for plane in &self.planes {
            match plane.ratio_to_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }
}

impl<T, const N: usize> Periodic<Cartessian<T, N>> for SimpleBox<T, N>
where
    T: Scalar,
    SimplePlanePair<T>: Periodic<Cartessian<T, N>>,
{
    fn find_pair(&self, pos: &Cartessian<T, N>) -> Cartessian<T, N> {
        let mut result = pos.clone();
        for planepair in &self.planes {
            planepair.find_pair_mut(&mut result);
        }
        return result;
    }

    fn find_pair_mut(&self, pos: &mut Cartessian<T, N>) {
        for planepair in &self.planes {
            planepair.find_pair_mut(pos);
        }
    }
}

impl<T, const N: usize> IntBoundary<Cartessian<T, N>> for SimpleBox<T, N>
where
    T: Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
    SimplePlanePair<T>: IntBoundary<Cartessian<T, N>>,
{
    fn check_inclusion(&self, pos: &Cartessian<T, N>) -> bool {
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos: &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos: &Cartessian<T, N>) -> Cartessian<T, N> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return v,
                None => {}
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<Cartessian<T, N>> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn find_intersect_unsafe(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<Cartessian<T, N>> {
        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }
}

impl<T, const N: usize> FloatBoundary<CartessianND<T>> for SimpleBox<T, N>
where
    T: Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
    SimplePlanePair<T>: FloatBoundary<CartessianND<T>>,
{
    fn check_inclusion(&self, pos: &CartessianND<T>) -> bool {
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos: &CartessianND<T>) -> Option<CartessianND<T>> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos: &CartessianND<T>) -> CartessianND<T> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return v,
                None => {}
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(
        &self,
        pos: &CartessianND<T>,
        movement: &CartessianND<T>,
    ) -> Option<CartessianND<T>> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn find_intersect_unsafe(
        &self,
        pos: &CartessianND<T>,
        movement: &CartessianND<T>,
    ) -> Option<CartessianND<T>> {
        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn ratio_to_intersect(&self, pos: &CartessianND<T>, movement: &CartessianND<T>) -> Option<T> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.ratio_to_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn ratio_to_intersect_unsafe(
        &self,
        pos: &CartessianND<T>,
        movement: &CartessianND<T>,
    ) -> Option<T> {
        for plane in &self.planes {
            match plane.ratio_to_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }
}

impl<T, const N: usize> Periodic<CartessianND<T>> for SimpleBox<T, N>
where
    T: Scalar,
    SimplePlanePair<T>: Periodic<CartessianND<T>>,
{
    fn find_pair(&self, pos: &CartessianND<T>) -> CartessianND<T> {
        let mut result = pos.clone();
        for planepair in &self.planes {
            planepair.find_pair_mut(&mut result);
        }
        return result;
    }

    fn find_pair_mut(&self, pos: &mut CartessianND<T>) {
        for planepair in &self.planes {
            planepair.find_pair_mut(pos);
        }
    }
}

impl<T, const N: usize> IntBoundary<CartessianND<T>> for SimpleBox<T, N>
where
    T: Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
    SimplePlanePair<T>: IntBoundary<CartessianND<T>>,
{
    fn check_inclusion(&self, pos: &CartessianND<T>) -> bool {
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos: &CartessianND<T>) -> Option<CartessianND<T>> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos: &CartessianND<T>) -> CartessianND<T> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return v,
                None => {}
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(
        &self,
        pos: &CartessianND<T>,
        movement: &CartessianND<T>,
    ) -> Option<CartessianND<T>> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn find_intersect_unsafe(
        &self,
        pos: &CartessianND<T>,
        movement: &CartessianND<T>,
    ) -> Option<CartessianND<T>> {
        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cube<V: Vector> {
    center: V,
    radius: <V as Vector>::Item,
}

impl<V: Vector> Cube<V> {
    pub fn new(center: V, radius: <V as Vector>::Item) -> Self
    where
        <V as Vector>::Item: Copy,
    {
        Self {
            center: center,
            radius: radius,
        }
    }
}

macro_rules! impl_clap_cube {
    (fixed, $ty : ident) => {
        impl<'h, const N : usize> CommandBuilder<'h, 2> for Cube<Cartessian<$ty, N>> {
            const SUBCOMMAND: &'h str = "Cube";

            fn args() -> [Arg<'h>; 2] {
                [Arg::new("center")
                    .short('c')
                    .long("center")
                    .value_name("CENTER")
                    .takes_value(true)
                    .help("Center of cube. format : 0,2,3"),
                 Arg::new("radius")
                    .short('r')
                    .long("radius")
                    .value_name("RADIUS")
                    .takes_value(true)
                    .value_parser(clap::value_parser!($ty))
                    .help("Radius of cube"),]
            }
        }

        impl<const N : usize> From<&ArgMatches> for Cube<Cartessian<$ty, N>> {
            fn from(m: &ArgMatches) -> Self {
                let center : Cartessian<$ty, N> = m.get_one::<String>("center").unwrap().parse().unwrap();
                let radius: $ty = *m.get_one::<$ty>("radius").unwrap();
                Cube::<Cartessian<$ty, N>>::new(center, radius)
            }
        }
    };
    (nd, $ty : ident) => {
        impl<'h> CommandBuilder<'h, 2> for Cube<CartessianND<$ty>> {
            const SUBCOMMAND: &'h str = "Cube";

            fn args() -> [Arg<'h>; 2] {
                [Arg::new("center")
                    .short('c')
                    .long("center")
                    .value_name("CENTER")
                    .takes_value(true)
                    .help("Center of cube. format : 0,2,3"),
                 Arg::new("radius")
                    .short('r')
                    .long("radius")
                    .value_name("RADIUS")
                    .takes_value(true)
                    .value_parser(clap::value_parser!($ty))
                    .help("Radius of cube"),]
            }
        }

        impl From<&ArgMatches> for Cube<CartessianND<$ty>> {
            fn from(m: &ArgMatches) -> Self {
                let center : CartessianND<$ty> = m.get_one::<String>("center").unwrap().parse().unwrap();
                let radius: $ty = *m.get_one::<$ty>("radius").unwrap();
                Cube::<CartessianND<$ty>>::new(center, radius)
            }
        }
    };
    ($ty : ident) => {
        impl_clap_cube!(fixed, $ty);
        impl_clap_cube!(nd, $ty);
    }
}

impl_clap_cube!(f32);
impl_clap_cube!(f64);
impl_clap_cube!(i8);
impl_clap_cube!(i16);
impl_clap_cube!(i32);
impl_clap_cube!(i64);
impl_clap_cube!(isize);



macro_rules! impl_float_cube {
    ($ty : ident) => {
        impl<const N: usize> FloatBoundary<Cartessian<$ty, N>> for Cube<Cartessian<$ty, N>> {
            fn check_inclusion(&self, pos: &Cartessian<$ty, N>) -> bool {
                for (c, x) in self.center.into_iter().zip(pos) {
                    if (*x - *c) > self.radius || (*x - *c) < -self.radius {
                        return false;
                    }
                }
                return true;
            }

            fn normal_at(&self, pos: &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let mut vec = Cartessian::<$ty, N>::default();
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon) {
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon) {
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let mut vec = Cartessian::<$ty, N>::default();
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon) {
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon) {
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                for idx in 0..N {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        return Some(pos + movement * t);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        return Some(pos + movement * t);
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                for idx in 0..N {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        return Some(pos + movement * t);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        return Some(pos + movement * t);
                    }
                }
                return None;
            }

            fn ratio_to_intersect(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<$ty> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                for idx in 0..N {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        return Some(t);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        return Some(t);
                    }
                }
                return None;
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<$ty> {
                for idx in 0..N {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        return Some(t);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        return Some(t);
                    }
                }
                return None;
            }
        }

        impl<const N: usize> Periodic<Cartessian<$ty, N>> for Cube<Cartessian<$ty, N>> {
            fn find_pair(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let mut result = pos.clone();
                result.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius;
                    }
                });
                return result;
            }

            fn find_pair_mut(&self, pos: &mut Cartessian<$ty, N>) {
                pos.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius;
                    }
                });
            }
        }

        impl FloatBoundary<CartessianND<$ty>> for Cube<CartessianND<$ty>> {
            fn check_inclusion(&self, pos: &CartessianND<$ty>) -> bool {
                if self.center.dim() != pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                for (c, x) in self.center.into_iter().zip(pos) {
                    if *x > *c + self.radius || *x < *c - self.radius {
                        return false;
                    }
                }
                return true;
            }

            fn normal_at(&self, pos: &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if self.center.dim() != pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon) {
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon) {
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                let epsilon = <$ty as AbsDiffEq>::default_epsilon();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if (*x - *c).abs_diff_eq(&self.radius, epsilon) {
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if (*x - *c).abs_diff_eq(&-self.radius, epsilon) {
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let n = pos.dim();
                for idx in 0..n {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                let n = pos.dim();
                for idx in 0..n {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        let mut result = pos.clone();
                        result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                        return Some(result);
                    }
                }
                return None;
            }

            fn ratio_to_intersect(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<$ty> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let n = pos.dim();
                for idx in 0..n {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        return Some(t);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        return Some(t);
                    }
                }
                return None;
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<$ty> {
                let n = pos.dim();
                for idx in 0..n {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius {
                        let t = (c - self.radius - p) / m;
                        return Some(t);
                    } else if p + m > c + self.radius {
                        let t = (c + self.radius - p) / m;
                        return Some(t);
                    }
                }
                return None;
            }
        }

        impl Periodic<CartessianND<$ty>> for Cube<CartessianND<$ty>> {
            fn find_pair(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let mut result = pos.clone();
                result.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius;
                    }
                });
                return result;
            }

            fn find_pair_mut(&self, pos: &mut CartessianND<$ty>) {
                pos.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius;
                    }
                });
            }
        }
    };
}

impl_float_cube!(f32);
impl_float_cube!(f64);

macro_rules! impl_int_cube {
    ($ty : ident) => {
        impl<const N: usize> IntBoundary<Cartessian<$ty, N>> for Cube<Cartessian<$ty, N>> {
            fn check_inclusion(&self, pos: &Cartessian<$ty, N>) -> bool {
                for (c, x) in self.center.into_iter().zip(pos) {
                    if (*x - *c) > self.radius || (*x - *c) < -self.radius {
                        return false;
                    }
                }
                return true;
            }

            fn normal_at(&self, pos: &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let mut vec = Cartessian::<$ty, N>::default();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if *x == *c + &self.radius {
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if *x == *c - &self.radius {
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let mut vec = Cartessian::<$ty, N>::default();
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if *x == *c + &self.radius {
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if *x == *c - &self.radius {
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                for idx in 0..N {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius || p + m > c + self.radius {
                        return Some(pos.clone());
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                for idx in 0..N {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius || p + m > c + self.radius {
                        return Some(pos.clone());
                    }
                }

                return None;
            }
        }

        impl<const N: usize> Periodic<Cartessian<$ty, N>> for Cube<Cartessian<$ty, N>> {
            fn find_pair(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let mut result = pos.clone();
                result.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius + 1 as $ty;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius + 1 as $ty;
                    }
                });
                return result;
            }

            fn find_pair_mut(&self, pos: &mut Cartessian<$ty, N>) {
                pos.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius + 1 as $ty;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius + 1 as $ty;
                    }
                });
            }
        }

        impl IntBoundary<CartessianND<$ty>> for Cube<CartessianND<$ty>> {
            fn check_inclusion(&self, pos: &CartessianND<$ty>) -> bool {
                if self.center.dim() != pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                for (c, x) in self.center.into_iter().zip(pos) {
                    if *x > *c + self.radius || *x < *c - self.radius {
                        return false;
                    }
                }
                return true;
            }

            fn normal_at(&self, pos: &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                if self.center.dim() != pos.dim() {
                    panic!("Dimensionality of plane and position vector are not compatible.");
                }

                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if *x == *c + &self.radius {
                        vec[idx] = -1 as $ty;
                        return Some(vec);
                    } else if *x == *c - &self.radius {
                        vec[idx] = 1 as $ty;
                        return Some(vec);
                    }
                }
                return None;
            }

            fn normal_at_unsafe(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let mut vec = CartessianND::<$ty>::zeros(pos.dim());
                for (idx, (c, x)) in self.center.into_iter().zip(pos).enumerate() {
                    if *x == *c + &self.radius {
                        vec[idx] = -1 as $ty;
                        return vec;
                    } else if *x == *c - &self.radius {
                        vec[idx] = 1 as $ty;
                        return vec;
                    }
                }
                panic!("Position is not on the boundaries");
            }

            fn find_intersect(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                if !self.check_inclusion(pos) {
                    panic!(
                        "State cannot live outside of system boundary. Move from outside occurs"
                    );
                }

                let n = pos.dim();
                for idx in 0..n {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius || p + m > c + self.radius {
                        return Some(pos.clone());
                    }
                }
                return None;
            }

            fn find_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                let n = pos.dim();
                for idx in 0..n {
                    let (p, m, c) = (pos[idx], movement[idx], self.center[idx]);
                    if p + m < c - self.radius || p + m > c + self.radius {
                        return Some(pos.clone());
                    }
                }
                return None;
            }
        }

        impl Periodic<CartessianND<$ty>> for Cube<CartessianND<$ty>> {
            fn find_pair(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let mut result = pos.clone();
                result.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius + 1 as $ty;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius + 1 as $ty;
                    }
                });
                return result;
            }

            fn find_pair_mut(&self, pos: &mut CartessianND<$ty>) {
                pos.zip_mut_with(&self.center, |x, y| {
                    if *x < *y - self.radius {
                        *x += 2 as $ty * self.radius + 1 as $ty;
                    } else if *y + self.radius < *x {
                        *x -= 2 as $ty * self.radius + 1 as $ty;
                    }
                });
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

#[derive(Clone, Debug, PartialEq)]
pub struct Parallelogram<T: Scalar, const N: usize> {
    pub planes: [PlanePair<Cartessian<T, N>>; N],
}

impl<T: Scalar, const N: usize> Parallelogram<T, N> {
    pub fn new(planes: [PlanePair<Cartessian<T, N>>; N]) -> Result<Self, Error>
    where
        T: Clone + AbsDiffEq,
        <T as AbsDiffEq>::Epsilon: Clone,
    {
        for i in 0..N - 1 {
            for j in i + 1..N {
                if planes[i]
                    .normal_vec
                    .abs_diff_eq(&planes[j].normal_vec, <T as AbsDiffEq>::default_epsilon())
                {
                    return Err(Error::make_error_syntax(
                        crate::prelude::ErrorCode::InvalidArgumentInput,
                    ));
                }
            }
        }

        Ok(Self {
            planes: planes.clone(),
        })
    }
}

impl<T, const N: usize> FloatBoundary<Cartessian<T, N>> for Parallelogram<T, N>
where
    T: Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq,
    PlanePair<Cartessian<T, N>>: FloatBoundary<Cartessian<T, N>>,
{
    fn check_inclusion(&self, pos: &Cartessian<T, N>) -> bool {
        self.planes.iter().all(|p| p.check_inclusion(pos))
    }

    fn normal_at(&self, pos: &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn normal_at_unsafe(&self, pos: &Cartessian<T, N>) -> Cartessian<T, N> {
        for plane in &self.planes {
            match plane.normal_at(pos) {
                Some(v) => return v,
                None => {}
            }
        }
        panic!("Position is not on the boundaries");
    }

    fn find_intersect(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<Cartessian<T, N>> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn find_intersect_unsafe(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<Cartessian<T, N>> {
        for plane in &self.planes {
            match plane.find_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn ratio_to_intersect(&self, pos: &Cartessian<T, N>, movement: &Cartessian<T, N>) -> Option<T> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        } else if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        for plane in &self.planes {
            match plane.ratio_to_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }

    fn ratio_to_intersect_unsafe(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<T> {
        for plane in &self.planes {
            match plane.ratio_to_intersect_unsafe(pos, movement) {
                Some(v) => return Some(v),
                None => {}
            }
        }
        return None;
    }
}

impl<T, const N: usize> Periodic<Cartessian<T, N>> for Parallelogram<T, N>
where
    T: Scalar,
    PlanePair<Cartessian<T, N>>: Periodic<Cartessian<T, N>>,
{
    fn find_pair(&self, pos: &Cartessian<T, N>) -> Cartessian<T, N> {
        let mut result = pos.clone();
        for planepair in &self.planes {
            planepair.find_pair_mut(&mut result);
        }
        return result;
    }

    fn find_pair_mut(&self, pos: &mut Cartessian<T, N>) {
        for planepair in &self.planes {
            planepair.find_pair_mut(pos);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::{Cartessian2D, Cartessian3D, CartessianND};
    use approx::assert_abs_diff_eq;
    use clap::Command;

    #[test]
    fn test_simpleplane() {
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
        assert_eq!(
            plane.find_intersect(&a, &movement),
            Some(Cartessian2D::new([1f64; 2]))
        );
        assert_eq!(
            plane.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([1f64; 2]))
        );

        movement[0] = -0.5f64;
        assert_eq!(plane.find_intersect(&a, &movement), None);
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), None);

        // float ND
        let mut a: CartessianND<f64> = CartessianND::new(vec![2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(
            plane.normal_at(&a),
            Some(CartessianND::new(vec![1f64, 0.0]))
        );
        assert_eq!(
            plane.normal_at_unsafe(&a),
            CartessianND::new(vec![1f64, 0.0])
        );

        a[0] = 2f64;
        let mut movement = CartessianND::new(vec![-2f64; 2]);
        assert_eq!(
            plane.find_intersect(&a, &movement),
            Some(CartessianND::new(vec![1f64; 2]))
        );
        assert_eq!(
            plane.find_intersect_unsafe(&a, &movement),
            Some(CartessianND::new(vec![1f64; 2]))
        );

        movement[0] = -0.5f64;
        assert_eq!(plane.find_intersect(&a, &movement), None);
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), None);

        let plane = SimplePlane::new(0, 1, Direction::Positive);
        // int 2D
        let mut a: Cartessian<i32, 2> = Cartessian2D::new([2; 2]);
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
        assert_eq!(
            plane.find_intersect(&a, &movement),
            Some(Cartessian2D::new([1, 2]))
        );
        assert_eq!(
            plane.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([1, 2]))
        );

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
        assert_eq!(
            plane.find_intersect(&a, &movement),
            Some(CartessianND::new(vec![1, 2]))
        );
        assert_eq!(
            plane.find_intersect_unsafe(&a, &movement),
            Some(CartessianND::new(vec![1, 2]))
        );

        movement[0] = 0;
        assert_eq!(plane.find_intersect(&a, &movement), None);
        assert_eq!(plane.find_intersect_unsafe(&a, &movement), None);
    }

    #[test]
    #[should_panic]
    fn test_simpleplane_panic() {
        let plane = SimplePlane::new(3, 1f64, Direction::Positive);

        let a = Cartessian2D::new([2f64; 2]);
        plane.check_inclusion(&a);
    }

    #[test]
    #[should_panic]
    fn test_simpleplane_find_itersect() {
        let plane = SimplePlane::new(0, 1f64, Direction::Positive);

        let a = Cartessian2D::new([0f64; 2]);
        let movement = Cartessian2D::new([1f64, 0f64]);
        plane.find_intersect(&a, &movement);
    }

    #[test]
    fn test_simpleplane_clap() {
        let arg = Command::new("test")
            .args(SimplePlane::<f64>::args())
            .get_matches_from(vec!["test", "--idx", "0", "--pos", "1", "--dir", "0"]);
        let simpleplane = SimplePlane::<f64>::from(&arg);
        assert_eq!(
            simpleplane,
            SimplePlane::<f64>::new(0, 1.0, Direction::Negative)
        )
    }

    #[test]
    fn test_plane() {
        // 2d
        let plane = Plane::new(Cartessian2D::new([1f64, 1f64]), 0f64);
        assert_abs_diff_eq!(plane.normal_vec, Cartessian2D::new([1.0 / 2f64.sqrt(); 2]));

        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -4f64;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -2f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(
            plane.normal_at(&a),
            Some(Cartessian2D::new([1.0 / 2f64.sqrt(); 2]))
        );
        assert_eq!(
            plane.normal_at_unsafe(&a),
            Cartessian2D::new([1.0 / 2f64.sqrt(); 2])
        );

        a[0] = 2f64;
        let movement = Cartessian2D::new([-3f64; 2]);
        assert_eq!(
            plane.find_intersect(&a, &movement),
            Some(Cartessian2D::new([0f64; 2]))
        );
        assert_eq!(
            plane.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([0f64; 2]))
        );

        // nd
        let plane = Plane::new(CartessianND::new(vec![1f64, 1f64]), 0f64);
        assert_abs_diff_eq!(
            plane.normal_vec,
            CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])
        );

        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -4f64;
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -2f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(
            plane.normal_at(&a),
            Some(CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]))
        );
        assert_eq!(
            plane.normal_at_unsafe(&a),
            CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])
        );

        a[0] = 2f64;
        let movement = CartessianND::new(vec![-3f64; 2]);
        assert_eq!(
            plane.find_intersect(&a, &movement),
            Some(CartessianND::new(vec![0f64; 2]))
        );
        assert_eq!(
            plane.find_intersect_unsafe(&a, &movement),
            Some(CartessianND::new(vec![0f64; 2]))
        );
    }

    #[test]
    #[should_panic]
    fn test_plane_panic() {
        let plane = Plane::new(CartessianND::new(vec![3f64; 3]), 0f64);

        let a = CartessianND::new(vec![2f64; 2]);
        plane.check_inclusion(&a);
    }

    #[test]
    #[should_panic]
    fn test_plane_intersect() {
        let plane = Plane::new(CartessianND::new(vec![3f64; 2]), 0f64);

        let a = CartessianND::new(vec![-2f64; 2]);
        let movement = CartessianND::new(vec![0f64; 2]);
        plane.find_intersect(&a, &movement);
    }

    #[test]
    fn test_plane_clap() {
        let arg = Command::new("test")
            .args(Plane::<Cartessian2D<f64>>::args())
            .get_matches_from(vec!["test", "--normal", "0,1", "--const", "0"]);
        let plane = Plane::<Cartessian2D<f64>>::from(&arg);
        assert_eq!(
            plane,
            Plane::<Cartessian2D<f64>>::new(Cartessian2D::<f64>::new([0.0, 1.0]), 0.0)
        )
    }

    #[test]
    fn test_simpleplanepair() {
        // Float
        let planepair = SimplePlanePair::new(0, [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(planepair.pos[0], 0f64);
        assert_abs_diff_eq!(planepair.pos[1], 1f64);

        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(Cartessian2D::new([-1f64, 0.0]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            Cartessian2D::new([-1f64, 0.0])
        );

        a[0] = 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(Cartessian2D::new([1f64, 0.0]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            Cartessian2D::new([1f64, 0.0])
        );

        a[0] = -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0.5f64;
        let mut movement = Cartessian2D::new([-2f64; 2]);
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(Cartessian2D::new([0f64, 1.5f64]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([0f64, 1.5f64]))
        );
        movement[0] = 2f64;
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(Cartessian2D::new([1f64, 1.5f64]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([1f64, 1.5f64]))
        );

        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(CartessianND::new(vec![-1f64, 0.0]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            CartessianND::new(vec![-1f64, 0.0])
        );

        a[0] = 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(CartessianND::new(vec![1f64, 0.0]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            CartessianND::new(vec![1f64, 0.0])
        );

        a[0] = -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0.5f64;
        let mut movement = CartessianND::new(vec![-2f64; 2]);
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(CartessianND::new(vec![0f64, 1.5f64]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(CartessianND::new(vec![0f64, 1.5f64]))
        );
        movement[0] = 2f64;
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(CartessianND::new(vec![1f64, 1.5f64]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(CartessianND::new(vec![1f64, 1.5f64]))
        );

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
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(a.clone())
        );

        a[0] = 2;
        movement[0] = 1;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(a.clone()));
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(a.clone())
        );

        let mut a = CartessianND::new(vec![3; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 2;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(CartessianND::new(vec![-1, 0]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            CartessianND::new(vec![-1, 0])
        );

        a[0] = 1;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![1, 0])));
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            CartessianND::new(vec![1, 0])
        );

        a[0] = -1;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0;
        let mut movement = CartessianND::new(vec![-1, 0]);
        assert_eq!(planepair.find_intersect(&a, &movement), Some(a.clone()));
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(a.clone())
        );

        a[0] = 2;
        movement[0] = 1;
        assert_eq!(planepair.find_intersect(&a, &movement), Some(a.clone()));
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(a.clone())
        );
    }

    #[test]
    fn test_simple_plane_pair_periodic() {
        // Float, fixed size
        let planepair = SimplePlanePair::new(0, [0f64, 4f64]).unwrap();

        let pos1 = Cartessian2D::new([-1f64, 2f64]);
        assert_eq!(planepair.find_pair(&pos1), Cartessian2D::new([3f64, 2f64]));

        let pos2 = Cartessian2D::new([4.5f64, 1f64]);
        assert_eq!(
            planepair.find_pair(&pos2),
            Cartessian2D::new([0.5f64, 1f64])
        );

        // Float, free size
        let planepair = SimplePlanePair::new(0, [0f64, 4f64]).unwrap();

        let pos1 = CartessianND::new(vec![-1f64, 2f64]);
        assert_eq!(
            planepair.find_pair(&pos1),
            CartessianND::new(vec![3f64, 2f64])
        );

        let pos2 = CartessianND::new(vec![4.5f64, 1f64]);
        assert_eq!(
            planepair.find_pair(&pos2),
            CartessianND::new(vec![0.5f64, 1f64])
        );

        // Integer, fixed size
        let planepair = SimplePlanePair::new(0, [0i32, 4i32]).unwrap();

        let pos1 = Cartessian2D::new([-1i32, 2i32]);
        assert_eq!(planepair.find_pair(&pos1), Cartessian2D::new([4i32, 2i32]));

        let pos2 = Cartessian2D::new([6i32, 1i32]);
        assert_eq!(planepair.find_pair(&pos2), Cartessian2D::new([1i32, 1i32]));

        // Integer, free size
        let planepair = SimplePlanePair::new(0, [0i32, 4i32]).unwrap();

        let pos1 = CartessianND::new(vec![-1i32, 2i32]);
        assert_eq!(
            planepair.find_pair(&pos1),
            CartessianND::new(vec![4i32, 2i32])
        );

        let pos2 = CartessianND::new(vec![6i32, 1i32]);
        assert_eq!(
            planepair.find_pair(&pos2),
            CartessianND::new(vec![1i32, 1i32])
        );
    }

    #[test]
    fn test_simpleplanepair_clap() {
        let arg = Command::new("test")
            .args(SimplePlanePair::<f64>::args())
            .get_matches_from(vec!["test", "--idx", "0", "--pos", "[1,2]"]);
        let simpleplanepair = SimplePlanePair::<f64>::from(&arg);
        assert_eq!(
            simpleplanepair,
            SimplePlanePair::<f64>::new(0, [1.0, 2.0]).unwrap()
        )
    }

    #[test]
    #[should_panic]
    fn test_planepair_panic() {
        let planepair = SimplePlanePair::new(0, [1f64, 0f64]).unwrap();
        let a = Cartessian2D::new([3f64; 2]);
        let movement = Cartessian2D::new([-2f64; 2]);
        planepair.find_intersect(&a, &movement);
    }

    #[test]
    fn test_planepair() {
        let planepair = PlanePair::new(Cartessian2D::new([1f64, 1f64]), [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(
            planepair.normal_vec,
            Cartessian2D::new([1.0 / 2f64.sqrt(); 2])
        );

        let mut a = Cartessian2D::new([1f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(Cartessian2D::new([-1.0 / 2f64.sqrt(); 2]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            Cartessian2D::new([-1.0 / 2f64.sqrt(); 2])
        );

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a *= -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(Cartessian2D::new([1.0 / 2f64.sqrt(); 2]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            Cartessian2D::new([1.0 / 2f64.sqrt(); 2])
        );

        let a = Cartessian2D::new([0f64, 0.5f64]);
        let mut movement = Cartessian2D::new([0f64, 1f64]);
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(Cartessian2D::new([0f64, 1f64]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([0f64, 1f64]))
        );
        movement[1] = -1f64;
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(Cartessian2D::new([0f64; 2]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([0f64; 2]))
        );

        let planepair = PlanePair::new(CartessianND::new(vec![1f64, 1f64]), [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(
            planepair.normal_vec,
            CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])
        );

        let mut a = CartessianND::new(vec![1f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(CartessianND::new(vec![-1.0 / 2f64.sqrt(); 2]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            CartessianND::new(vec![-1.0 / 2f64.sqrt(); 2])
        );

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a *= -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(
            planepair.normal_at(&a),
            Some(CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]))
        );
        assert_eq!(
            planepair.normal_at_unsafe(&a),
            CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])
        );

        let a = CartessianND::new(vec![0f64, 0.5f64]);
        let mut movement = CartessianND::new(vec![0f64, 1f64]);
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(CartessianND::new(vec![0f64, 1f64]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(CartessianND::new(vec![0f64, 1f64]))
        );
        movement[1] = -1f64;
        assert_eq!(
            planepair.find_intersect(&a, &movement),
            Some(CartessianND::new(vec![0f64; 2]))
        );
        assert_eq!(
            planepair.find_intersect_unsafe(&a, &movement),
            Some(CartessianND::new(vec![0f64; 2]))
        );
    }

    #[test]
    fn test_plane_pair_periodic() {
        // Float, fixed size
        let normal = Cartessian2D::new([1f64, 0f64]);
        let planepair = PlanePair::new(normal, [0f64, 4f64]).unwrap();

        let pos1 = Cartessian2D::new([-1f64, 2f64]);
        assert_eq!(planepair.find_pair(&pos1), Cartessian2D::new([3f64, 2f64]));

        let pos2 = Cartessian2D::new([4.5f64, 1f64]);
        assert_eq!(
            planepair.find_pair(&pos2),
            Cartessian2D::new([0.5f64, 1f64])
        );

        // Float, free size
        let normal = CartessianND::new(vec![1f64, 0f64]);
        let planepair = PlanePair::new(normal, [0f64, 4f64]).unwrap();

        let pos1 = CartessianND::new(vec![-1f64, 2f64]);
        assert_eq!(
            planepair.find_pair(&pos1),
            CartessianND::new(vec![3f64, 2f64])
        );

        let pos2 = CartessianND::new(vec![4.5f64, 1f64]);
        assert_eq!(
            planepair.find_pair(&pos2),
            CartessianND::new(vec![0.5f64, 1f64])
        );
    }

    #[test]
    fn test_planepair_clap() {
        let arg = Command::new("test")
            .args(PlanePair::<Cartessian2D<f64>>::args())
            .get_matches_from(vec!["test", "--normal", "0,1", "--const", "[1,2]"]);
        let simpleplanepair = PlanePair::<Cartessian2D<f64>>::from(&arg);
        assert_eq!(
            simpleplanepair,
            PlanePair::<Cartessian2D<f64>>::new(Cartessian2D::<f64>::new([0.0, 1.0]), [1.0, 2.0])
                .unwrap()
        )
    }

    #[test]
    fn test_simplebox() {
        let simplebox: SimpleBox<i32, 3> =
            SimpleBox::cube_with_center(&Cartessian3D::new([0; 3]), 3).unwrap();
        for plane in &simplebox.planes {
            match plane.idx {
                0 => {
                    assert_eq!(plane, &SimplePlanePair::new(0, [-3, 3]).unwrap());
                }
                1 => {
                    assert_eq!(plane, &SimplePlanePair::new(1, [-3, 3]).unwrap());
                }
                2 => {
                    assert_eq!(plane, &SimplePlanePair::new(2, [-3, 3]).unwrap());
                }
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

        a[1] = 0;
        vec[1] = 0;
        a[2] = 3;
        vec[2] = -1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a[2] = -3;
        vec[2] = 1;
        assert_eq!(simplebox.check_inclusion(&a), true);
        assert_eq!(simplebox.normal_at(&a), Some(vec.clone()));
        assert_eq!(simplebox.normal_at_unsafe(&a), vec.clone());

        a.coord = [3, 0, 0];
        vec.coord = [1, 0, 0];
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = -3;
        vec[0] = -1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = 0;
        vec[0] = 0;
        a[1] = 3;
        vec[1] = 1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = -3;
        vec[1] = -1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = 0;
        vec[1] = 0;
        a[2] = 3;
        vec[2] = 1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[2] = -3;
        vec[2] = -1;
        assert_eq!(simplebox.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(simplebox.find_intersect_unsafe(&a, &vec), Some(a.clone()));
    }

    #[test]
    fn test_simplebox_periodic() {
        // Float, Fixed size
        let simplebox: SimpleBox<f64, 3> =
            SimpleBox::cube_with_center(&Cartessian3D::new([0f64; 3]), 2.0).unwrap();

        let mut pos = Cartessian3D::new([0f64; 3]);
        let mut res = pos.clone();
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        pos[2] = -3f64;
        res[2] = 1f64;
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        // Float, Free size
        let simplebox: SimpleBox<f64, 3> =
            SimpleBox::cube_with_center(&CartessianND::new(vec![0f64; 3]), 2.0).unwrap();

        let mut pos = CartessianND::new(vec![0f64; 3]);
        let mut res = pos.clone();
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        pos[2] = -3f64;
        res[2] = 1f64;
        simplebox.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        // Int, Fixed size
        let simplebox: SimpleBox<i32, 3> =
            SimpleBox::cube_with_center(&Cartessian3D::new([0i32; 3]), 2).unwrap();

        let mut pos = Cartessian3D::new([0i32; 3]);
        let mut res = pos.clone();
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        pos[2] = -3i32;
        res[2] = 2i32;
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        // Int, Free size
        let simplebox: SimpleBox<i32, 3> =
            SimpleBox::cube_with_center(&CartessianND::new(vec![0i32; 3]), 2).unwrap();

        let mut pos = CartessianND::new(vec![0i32; 3]);
        let mut res = pos.clone();
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        pos[2] = -3i32;
        res[2] = 2i32;
        simplebox.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);
    }

    #[test]
    fn test_simplebox_clap() {
        let arg = Command::new("test")
            .args(SimpleBox::<f64, 2>::args())
            .get_matches_from(vec!["test", "--pairs", "[[1,2],[2,3]]"]);
        let simplebox = SimpleBox::<f64, 2>::from(&arg);
        assert_eq!(
            simplebox,
            SimpleBox::<f64, 2>::from_pairs([[1.0, 2.0], [2.0, 3.0]])
        )
    }

    #[test]
    fn test_cube() {
        let cube = Cube::new(Cartessian3D::new([0; 3]), 3);
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

        a[1] = 0;
        vec[1] = 0;
        a[2] = 3;
        vec[2] = -1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[2] = -3;
        vec[2] = 1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a.coord = [3, 0, 0];
        vec.coord = [1, 0, 0];
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = -3;
        vec[0] = -1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[0] = 0;
        vec[0] = 0;
        a[1] = 3;
        vec[1] = 1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = -3;
        vec[1] = -1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[1] = 0;
        vec[1] = 0;
        a[2] = 3;
        vec[2] = 1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));

        a[2] = -3;
        vec[2] = -1;
        assert_eq!(cube.find_intersect(&a, &vec), Some(a.clone()));
        assert_eq!(cube.find_intersect_unsafe(&a, &vec), Some(a.clone()));
    }

    #[test]
    fn test_cube_periodic() {
        // Float, Fixed size
        let cube: Cube<Cartessian3D<f64>> = Cube::new(Cartessian3D::new([0f64; 3]), 2.0);

        let mut pos = Cartessian3D::new([0f64; 3]);
        let mut res = pos.clone();
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        pos[2] = -3f64;
        res[2] = 1f64;
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        // Float, Free size
        let cube: Cube<CartessianND<f64>> = Cube::new(CartessianND::new(vec![0f64; 3]), 2.0);

        let mut pos = CartessianND::new(vec![0f64; 3]);
        let mut res = pos.clone();
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        pos[2] = -3f64;
        res[2] = 1f64;
        cube.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        // Int, Fixed size
        let cube: Cube<Cartessian3D<i32>> = Cube::new(Cartessian3D::new([0i32; 3]), 2);

        let mut pos = Cartessian3D::new([0i32; 3]);
        let mut res = pos.clone();
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        pos[2] = -3i32;
        res[2] = 2i32;
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        // Int, Free size
        let cube: Cube<CartessianND<i32>> = Cube::new(CartessianND::new(vec![0i32; 3]), 2);

        let mut pos = CartessianND::new(vec![0i32; 3]);
        let mut res = pos.clone();
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);

        pos[0] = -3i32;
        res[0] = 2i32;
        pos[1] = -3i32;
        res[1] = 2i32;
        pos[2] = -3i32;
        res[2] = 2i32;
        cube.find_pair_mut(&mut pos);
        assert_eq!(&pos, &res);
    }

    #[test]
    fn test_parallelogram() {
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
        assert_eq!(
            parallelogram.find_intersect(&a, &movement),
            Some(Cartessian2D::new([1.5f64, 1f64]))
        );
        assert_eq!(
            parallelogram.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([1.5f64, 1f64]))
        );

        movement.coord = [-1f64; 2];
        assert_eq!(
            parallelogram.find_intersect(&a, &movement),
            Some(Cartessian2D::new([0.5f64, 0f64]))
        );
        assert_eq!(
            parallelogram.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([0.5f64, 0f64]))
        );

        movement.coord = [2f64, 0f64];
        assert_eq!(
            parallelogram.find_intersect(&a, &movement),
            Some(Cartessian2D::new([1.5f64, 0.5f64]))
        );
        assert_eq!(
            parallelogram.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([1.5f64, 0.5f64]))
        );

        movement.coord = [-2f64, 0f64];
        assert_eq!(
            parallelogram.find_intersect(&a, &movement),
            Some(Cartessian2D::new([0.5f64, 0.5f64]))
        );
        assert_eq!(
            parallelogram.find_intersect_unsafe(&a, &movement),
            Some(Cartessian2D::new([0.5f64, 0.5f64]))
        );
    }

    #[test]
    fn test_parallelogram_periodic() {
        // Float, Fixed size
        let parallelogram: Parallelogram<f64, 3> = Parallelogram::new([
            PlanePair::new(Cartessian3D::new([1f64, 0f64, 0f64]), [-2f64, 2f64]).unwrap(),
            PlanePair::new(Cartessian3D::new([0f64, 1f64, 0f64]), [-2f64, 2f64]).unwrap(),
            PlanePair::new(Cartessian3D::new([0f64, 0f64, 1f64]), [-2f64, 2f64]).unwrap(),
        ])
        .unwrap();

        let mut pos = Cartessian3D::new([0f64; 3]);
        let mut res = pos.clone();
        parallelogram.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        parallelogram.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        parallelogram.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);

        pos[0] = -3f64;
        res[0] = 1f64;
        pos[1] = -3f64;
        res[1] = 1f64;
        pos[2] = -3f64;
        res[2] = 1f64;
        parallelogram.find_pair_mut(&mut pos);
        assert_abs_diff_eq!(&pos, &res, epsilon = 1e-3);
    }

    #[test]
    fn test_serde_planes() {
        use serde_json::{from_str, to_string};

        // Direction
        let dir = Direction::Positive;
        let expected = r#""Positive""#;
        assert_eq!(expected, to_string(&dir).unwrap());

        let expected = from_str(&expected).unwrap();
        assert_eq!(dir, expected);

        // SimplePlane
        let simpleplane = SimplePlane::<f64>::new(0, 1.0, Direction::Positive);
        let expected = r#"{"idx":0,"pos":1.0,"dir_in":"Positive"}"#;
        assert_eq!(expected, to_string(&simpleplane).unwrap());

        let expected = from_str(&expected).unwrap();
        assert_eq!(simpleplane, expected);

        // Plane
        let normal_vec = Cartessian2D::new([0.0, 1.0]);
        let constant = 1.0f64;
        let plane: Plane<Cartessian2D<f64>> = Plane::new(normal_vec, constant);
        let expected = r#"{"normal_vec":{"coord":[0.0,1.0]},"constant":1.0}"#;
        assert_eq!(expected, to_string(&plane).unwrap());

        let expected = from_str(&expected).unwrap();
        assert_eq!(plane, expected);

        // Simple Plane Pair
        let simpleplanepair = SimplePlanePair::<f64>::new(0, [0.0, 1.0]).unwrap();
        let expected = r#"{"idx":0,"pos":[0.0,1.0]}"#;
        assert_eq!(expected, to_string(&simpleplanepair).unwrap());

        let expected = from_str(&expected).unwrap();
        assert_eq!(simpleplanepair, expected);

        // Plane Pair
        let normal_vec = Cartessian2D::new([0.0, 1.0]);
        let constant = [0.0f64, 1.0f64];
        let planepair: PlanePair<Cartessian2D<f64>> = PlanePair::new(normal_vec, constant).unwrap();
        let expected = r#"{"normal_vec":{"coord":[0.0,1.0]},"constant":[0.0,1.0]}"#;
        assert_eq!(expected, to_string(&planepair).unwrap());

        let expected = from_str(&expected).unwrap();
        assert_eq!(planepair, expected);

        // Cube
        let center = Cartessian2D::new([0.0, 0.0]);
        let radius = 1.0;
        let cube: Cube<Cartessian2D<f64>> = Cube::new(&center, radius);
        let expected = r#"{"center":{"coord":[0.0,0.0]},"radius":1.0}"#;
        assert_eq!(expected, to_string(&cube).unwrap());

        let expected = from_str(&expected).unwrap();
        assert_eq!(cube, expected);
    }
}
