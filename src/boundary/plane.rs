

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
pub struct SimplePlane<T>{
    idx : usize,
    pos : T,
    dir_out : bool,
}

impl<T> SimplePlane<T>{
    pub fn new(idx : usize, pos : T, dir_out : bool) -> Self{
        Self{
            idx,
            pos,
            dir_out,
        }
    }
}

impl<T, const N : usize> Boundary<Cartessian<T, N>> for SimplePlane<T>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &Cartessian<T, N>) -> bool{
        if self.idx >= N{
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        if self.dir_out{
            (self.pos - pos[self.idx]) >= T::zero()
        } else {
            (self.pos - pos[self.idx]) <= T::zero()
        }
    }

    fn normal_at(&self, pos : &Cartessian<T, N>) -> Option<Cartessian<T, N>>{
        if self.idx >= N{
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        if self.pos.abs_diff_eq(&pos[self.idx], <T as AbsDiffEq>::default_epsilon()){
            let mut normal_vec = Cartessian::<T, N>::zeros();
            if self.dir_out{
                normal_vec[self.idx] = T::one();
            } else {
                normal_vec[self.idx] = -(T::one());
            }
            return Some(normal_vec)
        } else {
            return None;
        }
    }
}

impl<T> Boundary<CartessianND<T>> for SimplePlane<T>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &CartessianND<T>) -> bool{
        if self.idx >= pos.dim(){
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        if self.dir_out{
            (self.pos - pos[self.idx]) >= T::zero()
        } else {
            (self.pos - pos[self.idx]) <= T::zero()
        }
    }

    fn normal_at(&self, pos : &CartessianND<T>) -> Option<CartessianND<T>>{
        if self.idx >= pos.dim(){
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        if self.pos.abs_diff_eq(&pos[self.idx], <T as AbsDiffEq>::default_epsilon()){
            let mut normal_vec = CartessianND::zeros(pos.dim());
            if self.dir_out{
                normal_vec[self.idx] = T::one();
            } else {
                normal_vec[self.idx] = -(T::one());
            }
            return Some(normal_vec)
        } else {
            return None;
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Plane<V : Vector>{
    normal_vec : V,
    constant : <V as Vector>::Item,
}

impl<'a, V> Plane<V>
    where V : Vector + Clone + Norm + Div<<V as Norm>::Output, Output = V> + Dot<&'a V, Output = <V as Vector>::Item> + 'a,
          <V as Vector>::Item : Clone{
    pub fn new(normal_vec : &V, constant : <V as Vector>::Item) -> Self{
        let n = normal_vec.norm_l2();
        Self{
            normal_vec : (*normal_vec).clone() / n,
            constant,
        }
    }

    pub fn from_point(normal_vec : &'a V, point : &'a V) -> Self{
        let n = normal_vec.norm_l2();
        let vec = normal_vec.clone() / n;
        let c = vec.dot(point);
        Self{
            normal_vec : vec,
            constant : c,
        }
    }
}

impl<T, const N : usize> Boundary<Cartessian<T, N>> for Plane<Cartessian<T, N>>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &Cartessian<T, N>) -> bool{
        let d = self.normal_vec.dot(pos);

        d <= self.constant
    }

    fn normal_at(&self, pos : &Cartessian<T, N>) -> Option<Cartessian<T, N>>{
        let d = self.normal_vec.dot(pos);
        if self.constant.abs_diff_eq(&d, <T as AbsDiffEq>::default_epsilon()){
            return Some(self.normal_vec.clone());
        } else {
            return None;
        }
    }
}

impl<T> Boundary<CartessianND<T>> for Plane<CartessianND<T>>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &CartessianND<T>) -> bool{
        let d = self.normal_vec.dot(pos);

        d <= self.constant
    }

    fn normal_at(&self, pos : &CartessianND<T>) -> Option<CartessianND<T>>{
        let d = self.normal_vec.dot(pos);
        if self.constant.abs_diff_eq(&d, <T as AbsDiffEq>::default_epsilon()){
            return Some(self.normal_vec.clone());
        } else {
            return None;
        }
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

impl<T, const N : usize> Boundary<Cartessian<T, N>> for SimplePlanePair<T>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &Cartessian<T, N>) -> bool{
        if self.idx >= N{
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        self.pos[0] <= pos[self.idx]  && pos[self.idx] <= self.pos[1]
    }

    fn normal_at(&self, pos : &Cartessian<T, N>) -> Option<Cartessian<T, N>>{
        if self.idx >= N{
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        let mut normal_vec = Cartessian::<T, N>::zeros();
        if self.pos[0].abs_diff_eq(&pos[self.idx], <T as AbsDiffEq>::default_epsilon()){
            normal_vec[self.idx] = -T::one();
            return Some(normal_vec);
        } else if self.pos[1].abs_diff_eq(&pos[self.idx], <T as AbsDiffEq>::default_epsilon()){
            normal_vec[self.idx] = T::one();
            return Some(normal_vec);
        } else {
            return None;
        }
    }
}

impl<T> Boundary<CartessianND<T>> for SimplePlanePair<T>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &CartessianND<T>) -> bool{
        if self.idx >= pos.dim(){
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        self.pos[0] <= pos[self.idx]  && pos[self.idx] <= self.pos[1]
    }

    fn normal_at(&self, pos : &CartessianND<T>) -> Option<CartessianND<T>>{
        if self.idx >= pos.dim(){
            panic!("Dimensionality of plane and position vector are not compatible. index out of bounds");
        }

        let mut normal_vec = CartessianND::zeros(pos.dim());
        if self.pos[0].abs_diff_eq(&pos[self.idx], <T as AbsDiffEq>::default_epsilon()){
            normal_vec[self.idx] = -T::one();
            return Some(normal_vec);
        } else if self.pos[1].abs_diff_eq(&pos[self.idx], <T as AbsDiffEq>::default_epsilon()){
            normal_vec[self.idx] = T::one();
            return Some(normal_vec);
        } else {
            return None;
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PlanePair<V : Vector>{
    normal_vec : V,
    constant : [<V as Vector>::Item; 2],
}

impl<V> PlanePair<V>
    where V : Vector + Clone + Norm + Div<<V as Norm>::Output, Output = V>,
         <V as Vector>::Item : PartialOrd + AbsDiffEq{
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
            constant : [a, b],
        })
    }
}

impl<T, const N : usize> Boundary<Cartessian<T, N>> for PlanePair<Cartessian<T, N>>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &Cartessian<T, N>) -> bool{
        let d = self.normal_vec.dot(pos);

        self.constant[0] <= d && d <= self.constant[1]
    }

    fn normal_at(&self, pos : &Cartessian<T, N>) -> Option<Cartessian<T, N>>{
        let d = self.normal_vec.dot(pos);

        if self.constant[0].abs_diff_eq(&d, <T as AbsDiffEq>::default_epsilon()){
            return Some(-self.normal_vec.clone());
        } else if self.constant[1].abs_diff_eq(&d, <T as AbsDiffEq>::default_epsilon()){
            return Some(self.normal_vec.clone());
        } else {
            return None;
        }
    }
}

impl<T> Boundary<CartessianND<T>> for PlanePair<CartessianND<T>>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &CartessianND<T>) -> bool{
        let d = self.normal_vec.dot(pos);

        self.constant[0] <= d && d <= self.constant[1]
    }

    fn normal_at(&self, pos : &CartessianND<T>) -> Option<CartessianND<T>>{
        let d = self.normal_vec.dot(pos);

        if self.constant[0].abs_diff_eq(&d, <T as AbsDiffEq>::default_epsilon()){
            return Some(-self.normal_vec.clone());
        } else if self.constant[1].abs_diff_eq(&d, <T as AbsDiffEq>::default_epsilon()){
            return Some(self.normal_vec.clone());
        } else {
            return None;
        }
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

        if center.dim() != N {
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
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &Cartessian<T, N>) -> bool{
        self.planes.iter().map(|p| p.check_inclusion(pos)).all(|x| x)
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
}

impl<T, const N : usize> Boundary<CartessianND<T>> for SimpleBox<T, N>
    where T : Scalar + Neg<Output = T> + PartialOrd + AbsDiffEq{
    fn check_inclusion(&self, pos : &CartessianND<T>) -> bool{
        self.planes.iter().map(|p| p.check_inclusion(pos)).all(|x| x)
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
}

// #[derive(Clone, Debug, PartialEq)]
// pub struct Parallelogram<T, const N : usize>{
//     planes : [PlanePair<V>; N],
// }



#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::{Cartessian2D, Cartessian3D, CartessianND};
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_simpleplane(){
        let plane = SimplePlane::new(0, 1f64, true);

        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(Cartessian2D::new([1f64, 0.0])));

        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(CartessianND::new(vec![1f64, 0.0])));
    }

    #[test]
    #[should_panic]
    fn test_simpleplane_panic(){
        let plane = SimplePlane::new(3, 1f64, true);

        let a = Cartessian2D::new([2f64; 2]);
        plane.check_inclusion(&a);
    }

    #[test]
    fn test_plane(){
        let plane = Plane::new(&Cartessian2D::new([1f64, 1f64]), 0f64);
        assert_abs_diff_eq!(plane.normal_vec, Cartessian2D::new([1.0 / 2f64.sqrt(); 2]));

        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -4f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -2f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(Cartessian2D::new([1.0 / 2f64.sqrt(); 2])));

        let plane = Plane::new(&CartessianND::new(vec![1f64, 1f64]), 0f64);
        assert_abs_diff_eq!(plane.normal_vec, CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]));

        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(plane.check_inclusion(&a), false);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -4f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), None);

        a[0] = -2f64;
        assert_eq!(plane.check_inclusion(&a), true);
        assert_eq!(plane.normal_at(&a), Some(CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])));
    }

    #[test]
    #[should_panic]
    fn test_plane_panic(){
        let plane = Plane::new(&CartessianND::new(vec![3f64; 3]), 0f64);

        let a = CartessianND::new(vec![2f64; 2]);
        plane.check_inclusion(&a);
    }

    #[test]
    fn test_simpleplanepair(){
        let planepair = SimplePlanePair::new(0, [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(planepair.pos[0], 0f64);
        assert_abs_diff_eq!(planepair.pos[1], 1f64);

        let mut a = Cartessian2D::new([2f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([1f64, 0.0])));

        a[0] = 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([-1f64, 0.0])));

        a[0] = -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);


        let mut a = CartessianND::new(vec![2f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![1f64, 0.0])));

        a[0] = 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a[0] = 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![-1f64, 0.0])));

        a[0] = -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);
    }

    #[test]
    fn test_planepair(){
        let planepair = PlanePair::new(Cartessian2D::new([1f64, 1f64]), [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(planepair.normal_vec, Cartessian2D::new([1.0 / 2f64.sqrt(); 2]));


        let mut a = Cartessian2D::new([1f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 1.0/2.0f64.sqrt();
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(Cartessian2D::new([1.0 / 2f64.sqrt(); 2])));

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a *= -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(-Cartessian2D::new([1.0 / 2f64.sqrt(); 2])));

        let planepair = PlanePair::new(CartessianND::new(vec![1f64, 1f64]), [1f64, 0f64]).unwrap();
        assert_abs_diff_eq!(planepair.normal_vec, CartessianND::new(vec![1.0 / 2f64.sqrt(); 2]));

        let mut a = CartessianND::new(vec![1f64; 2]);
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 1.0/2.0f64.sqrt();
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])));

        a *= 0.5f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), None);

        a *= -1f64;
        assert_eq!(planepair.check_inclusion(&a), false);
        assert_eq!(planepair.normal_at(&a), None);

        a *= 0f64;
        assert_eq!(planepair.check_inclusion(&a), true);
        assert_eq!(planepair.normal_at(&a), Some(-CartessianND::new(vec![1.0 / 2f64.sqrt(); 2])));
    }

    #[test]
    fn test_cube(){
        let cube = SimpleBox::cube_with_center(&Cartessian3D::new([0; 3]), 3).unwrap();
        for plane in &cube.planes{
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

        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), None);

        a[0] = -3;
        vec[0] = -1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[0] = 3;
        vec[0] = 1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[0] = 5;
        assert_eq!(cube.check_inclusion(&a), false);
        assert_eq!(cube.normal_at(&a), None);

        a[0] = 0;
        a[1] = 3;
        vec[0] = 0;
        vec[1] = 1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[1] = -3;
        vec[1] = -1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[1] = 0; vec[1] = 0;
        a[2] = 3; vec[2] = 1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));

        a[2] = -3; vec[2] = -1;
        assert_eq!(cube.check_inclusion(&a), true);
        assert_eq!(cube.normal_at(&a), Some(vec.clone()));
    }

}
