use super::{FloatBoundary, IntBoundary, NonPeriodic};
use crate::boundary::Periodic;
use crate::vector::basic::Map;
use crate::vector::product::{InnerProduct, Norm};
use crate::vector::Scalar;
use crate::vector::{Cartessian, CartessianND, Vector};
use approx::AbsDiffEq;
use serde::{Deserialize, Serialize};
use std::ops::Neg;

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Sphere<V: Vector> {
    center: V,
    radius: <V as Vector>::Item,
}

impl<V: Vector> Sphere<V>
where
    V: Clone,
    <V as Vector>::Item: Clone,
{
    pub fn new(center: &V, radius: <V as Vector>::Item) -> Self {
        Self {
            center: center.clone(),
            radius: radius.clone(),
        }
    }
}

macro_rules! impl_float_sphere {
    ($ty : ident) => {
        impl<const N: usize> FloatBoundary<Cartessian<$ty, N>> for Sphere<Cartessian<$ty, N>> {
            fn check_inclusion(&self, pos: &Cartessian<$ty, N>) -> bool {
                (pos - &self.center).norm_l2() <= self.radius
            }

            fn normal_at(&self, pos: &Cartessian<$ty, N>) -> Option<Cartessian<$ty, N>> {
                let vec_dr = &self.center - pos;
                let dr = vec_dr.norm_l2();

                if dr.abs_diff_eq(&self.radius, <$ty as AbsDiffEq>::default_epsilon()) {
                    return Some(vec_dr / dr);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let vec_dr = &self.center - pos;
                let dr = vec_dr.norm_l2();
                vec_dr / dr
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
                } else if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let from_center = pos - &self.center;
                let (x2, dx2, xdx, r) = (
                    from_center.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    from_center.inner_product(movement),
                    self.radius,
                );
                let t = (-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2;
                return Some(pos + movement * t);
            }

            fn find_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<Cartessian<$ty, N>> {
                if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let from_center = pos - &self.center;
                let (x2, dx2, xdx, r) = (
                    from_center.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    from_center.inner_product(movement),
                    self.radius,
                );
                let t = (-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2;
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
                } else if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let from_center = pos - &self.center;
                let (x2, dx2, xdx, r) = (
                    from_center.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    from_center.inner_product(movement),
                    self.radius,
                );
                Some((-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2)
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &Cartessian<$ty, N>,
                movement: &Cartessian<$ty, N>,
            ) -> Option<$ty> {
                if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let from_center = pos - &self.center;
                let (x2, dx2, xdx, r) = (
                    from_center.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    from_center.inner_product(movement),
                    self.radius,
                );
                Some((-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2)
            }
        }

        impl<const N: usize> Periodic<Cartessian<$ty, N>> for Sphere<Cartessian<$ty, N>> {
            fn find_pair(&self, pos: &Cartessian<$ty, N>) -> Cartessian<$ty, N> {
                let r = pos.norm_l2();
                if r > self.radius {
                    return ((r - 2 as $ty * self.radius) / r) * pos;
                } else {
                    return pos.clone();
                }
            }

            fn find_pair_mut(&self, pos: &mut Cartessian<$ty, N>) {
                let r = pos.norm_l2();
                if r > self.radius {
                    let c = (r - 2 as $ty * self.radius) / r;
                    pos.map_inplace(|x| *x = *x * c);
                }
            }
        }

        impl FloatBoundary<CartessianND<$ty>> for Sphere<CartessianND<$ty>> {
            fn check_inclusion(&self, pos: &CartessianND<$ty>) -> bool {
                (pos - &self.center).norm_l2() <= self.radius
            }

            fn normal_at(&self, pos: &CartessianND<$ty>) -> Option<CartessianND<$ty>> {
                let vec_dr = &self.center - pos;
                let dr = vec_dr.norm_l2();

                if dr.abs_diff_eq(&self.radius, <$ty as AbsDiffEq>::default_epsilon()) {
                    return Some(vec_dr / dr);
                } else {
                    return None;
                }
            }

            fn normal_at_unsafe(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let vec_dr = &self.center - pos;
                let dr = vec_dr.norm_l2();
                vec_dr / dr
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
                } else if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let mut result = pos.clone();
                result.zip_mut_with(&self.center, |x, y| *x = *x - *y);
                let (x2, dx2, xdx, r) = (
                    result.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    result.inner_product(movement),
                    self.radius,
                );
                let t = (-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2;
                result.clone_from(pos);
                result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
                return Some(result);
            }

            fn find_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<CartessianND<$ty>> {
                if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let mut result = pos.clone();
                result.zip_mut_with(&self.center, |x, y| *x = *x - *y);
                let (x2, dx2, xdx, r) = (
                    result.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    result.inner_product(movement),
                    self.radius,
                );
                let t = (-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2;
                result.clone_from(pos);
                result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
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
                } else if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let from_center = pos - &self.center;
                let (x2, dx2, xdx, r) = (
                    from_center.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    from_center.inner_product(movement),
                    self.radius,
                );
                Some((-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2)
            }

            fn ratio_to_intersect_unsafe(
                &self,
                pos: &CartessianND<$ty>,
                movement: &CartessianND<$ty>,
            ) -> Option<$ty> {
                if self.check_inclusion(&(pos + movement)) {
                    return None;
                }

                let from_center = pos - &self.center;
                let (x2, dx2, xdx, r) = (
                    from_center.norm_l2_sqr(),
                    movement.norm_l2_sqr(),
                    from_center.inner_product(movement),
                    self.radius,
                );
                Some((-xdx + (xdx * xdx + dx2 * (r * r - x2)).sqrt()) / dx2)
            }
        }

        impl Periodic<CartessianND<$ty>> for Sphere<CartessianND<$ty>> {
            fn find_pair(&self, pos: &CartessianND<$ty>) -> CartessianND<$ty> {
                let r = pos.norm_l2();
                if r > self.radius {
                    return ((r - 2 as $ty * self.radius) / r) * pos;
                } else {
                    return pos.clone();
                }
            }

            fn find_pair_mut(&self, pos: &mut CartessianND<$ty>) {
                let r = pos.norm_l2();
                if r > self.radius {
                    let c = (r - 2 as $ty * self.radius) / r;
                    pos.map_inplace(|x| *x = *x * c);
                }
            }
        }
    };
}

impl_float_sphere!(f32);
impl_float_sphere!(f64);

#[derive(Clone, Debug)]
pub struct TaxiSphere<V: Vector> {
    center: V,
    radius: <V as Vector>::Item,
}

impl<V: Vector> TaxiSphere<V>
where
    V: Clone,
    <V as Vector>::Item: Clone,
{
    pub fn new(center: &V, radius: <V as Vector>::Item) -> Self {
        Self {
            center: center.clone(),
            radius: radius.clone(),
        }
    }
}

impl<V: Vector> NonPeriodic for TaxiSphere<V> {}

impl<T, const N: usize> IntBoundary<Cartessian<T, N>> for TaxiSphere<Cartessian<T, N>>
where
    T: Scalar + AbsDiffEq + PartialOrd + Neg<Output = T>,
    Cartessian<T, N>: Vector<Item = T> + Norm<Output = T>,
{
    fn check_inclusion(&self, pos: &Cartessian<T, N>) -> bool {
        (pos - &self.center).norm_l1() <= self.radius
    }

    fn normal_at(&self, pos: &Cartessian<T, N>) -> Option<Cartessian<T, N>> {
        let vec_dr = &self.center - pos;
        let dr = vec_dr.norm_l1();

        if dr.abs_diff_eq(&self.radius, <T as AbsDiffEq>::default_epsilon()) {
            let mut normal_vec = Cartessian::<T, N>::zeros();
            for (v, nv) in vec_dr.iter().zip(normal_vec.iter_mut()) {
                if *v > T::zero() {
                    *nv = T::one();
                } else if *v < T::zero() {
                    *nv = -T::one();
                }
            }
            let n = normal_vec.norm_l2();

            return Some(normal_vec / n);
        } else {
            return None;
        }
    }

    fn normal_at_unsafe(&self, pos: &Cartessian<T, N>) -> Cartessian<T, N> {
        let vec_dr = &self.center - pos;
        let dr = vec_dr.norm_l1();

        if dr.abs_diff_eq(&self.radius, <T as AbsDiffEq>::default_epsilon()) {
            let mut normal_vec = Cartessian::<T, N>::zeros();
            for (v, nv) in vec_dr.iter().zip(normal_vec.iter_mut()) {
                if *v > T::zero() {
                    *nv = T::one();
                } else if *v < T::zero() {
                    *nv = -T::one();
                }
            }
            let n = normal_vec.norm_l2();

            return normal_vec / n;
        }
        panic!("Position is not on the boundary");
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

        let (mut sum_dx, mut sum_x) = (T::zero(), self.radius);
        for idx in 0..N {
            let (x, dx, c) = (pos[idx], movement[idx], self.center[idx]);
            if x + dx > c {
                sum_dx += dx;
                sum_x -= x + c;
            } else {
                sum_dx -= dx;
                sum_x += x + c;
            }
        }
        let t = sum_x / sum_dx;
        return Some(pos + movement * t);
    }

    fn find_intersect_unsafe(
        &self,
        pos: &Cartessian<T, N>,
        movement: &Cartessian<T, N>,
    ) -> Option<Cartessian<T, N>> {
        if self.check_inclusion(&(pos + movement)) {
            return None;
        }

        let (mut sum_dx, mut sum_x) = (T::zero(), self.radius);
        for idx in 0..N {
            let (x, dx, c) = (pos[idx], movement[idx], self.center[idx]);
            if x + dx > c {
                sum_dx += dx;
                sum_x -= x + c;
            } else {
                sum_dx -= dx;
                sum_x += x + c;
            }
        }
        let t = sum_x / sum_dx;
        return Some(pos + movement * t);
    }
}

impl<T> IntBoundary<CartessianND<T>> for TaxiSphere<CartessianND<T>>
where
    T: Scalar + AbsDiffEq + PartialOrd + Neg<Output = T>,
    CartessianND<T>: Vector<Item = T> + Norm<Output = T>,
{
    fn check_inclusion(&self, pos: &CartessianND<T>) -> bool {
        (pos - &self.center).norm_l1() <= self.radius
    }

    fn normal_at(&self, pos: &CartessianND<T>) -> Option<CartessianND<T>> {
        let vec_dr = &self.center - pos;
        let dr = vec_dr.norm_l1();

        if dr.abs_diff_eq(&self.radius, <T as AbsDiffEq>::default_epsilon()) {
            let mut normal_vec = CartessianND::<T>::zeros(pos.dim());
            for (v, nv) in vec_dr.iter().zip(normal_vec.iter_mut()) {
                if *v > T::zero() {
                    *nv = T::one();
                } else if *v < T::zero() {
                    *nv = -T::one();
                }
            }
            let n = normal_vec.norm_l2();

            return Some(normal_vec / n);
        } else {
            return None;
        }
    }

    fn normal_at_unsafe(&self, pos: &CartessianND<T>) -> CartessianND<T> {
        let vec_dr = &self.center - pos;
        let dr = vec_dr.norm_l1();

        if dr.abs_diff_eq(&self.radius, <T as AbsDiffEq>::default_epsilon()) {
            let mut normal_vec = CartessianND::<T>::zeros(pos.dim());
            for (v, nv) in vec_dr.iter().zip(normal_vec.iter_mut()) {
                if *v > T::zero() {
                    *nv = T::one();
                } else if *v < T::zero() {
                    *nv = -T::one();
                }
            }
            let n = normal_vec.norm_l2();

            return normal_vec / n;
        }
        panic!("Position is not on the boundary");
    }

    fn find_intersect(
        &self,
        pos: &CartessianND<T>,
        movement: &CartessianND<T>,
    ) -> Option<CartessianND<T>> {
        if !self.check_inclusion(pos) {
            panic!("State cannot live outside of system boundary. Move from outside occurs");
        }

        let mut result = pos.clone();
        result.zip_mut_with(movement, |x, y| *x = *x + *y);
        if self.check_inclusion(&result) {
            return None;
        }

        let n = pos.dim();
        let (mut sum_dx, mut sum_x) = (T::zero(), self.radius);
        for idx in 0..n {
            let (x, dx, c) = (pos[idx], movement[idx], self.center[idx]);
            if x + dx > c {
                sum_dx += dx;
                sum_x -= x + c;
            } else {
                sum_dx -= dx;
                sum_x += x + c;
            }
        }
        let t = sum_x / sum_dx;
        result.clone_from(pos);
        result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
        return Some(result);
    }

    fn find_intersect_unsafe(
        &self,
        pos: &CartessianND<T>,
        movement: &CartessianND<T>,
    ) -> Option<CartessianND<T>> {
        let mut result = pos.clone();
        result.zip_mut_with(movement, |x, y| *x = *x + *y);
        if self.check_inclusion(&result) {
            return None;
        }

        let n = pos.dim();
        let (mut sum_dx, mut sum_x) = (T::zero(), self.radius);
        for idx in 0..n {
            let (x, dx, c) = (pos[idx], movement[idx], self.center[idx]);
            if x + dx > c {
                sum_dx += dx;
                sum_x -= x + c;
            } else {
                sum_dx -= dx;
                sum_x += x + c;
            }
        }
        let t = sum_x / sum_dx;
        result.clone_from(pos);
        result.zip_mut_with(movement, |x, y| *x = *x + *y * t);
        return Some(result);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::{Cartessian1D, Cartessian2D, Cartessian3D, CartessianND};
    use approx::assert_abs_diff_eq;
    use serde_json::{from_str, to_string};

    #[test]
    fn test_sphere1d() {
        let range = Sphere::new(&Cartessian1D::new([0f64]), 1f64);

        let mut a = Cartessian1D::new([-2f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), Cartessian1D::new([1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), Cartessian1D::new([1f64]));

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), Cartessian1D::new([-1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), Cartessian1D::new([-1f64]));

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 0f64;
        let mut movement = Cartessian1D::new([2f64]);
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            Cartessian1D::new([1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            Cartessian1D::new([1f64])
        );

        movement[0] = -2f64;
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            Cartessian1D::new([-1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            Cartessian1D::new([-1f64])
        );

        let range = Sphere::new(&CartessianND::new(vec![0f64]), 1f64);

        let mut a = CartessianND::new(vec![-2f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), CartessianND::new(vec![1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), CartessianND::new(vec![1f64]));

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), CartessianND::new(vec![-1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), CartessianND::new(vec![-1f64]));

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 0f64;
        let mut movement = CartessianND::new(vec![2f64]);
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            CartessianND::new(vec![1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            CartessianND::new(vec![1f64])
        );

        movement[0] = -2f64;
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            CartessianND::new(vec![-1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            CartessianND::new(vec![-1f64])
        );
    }

    #[test]
    fn test_sphere2d() {
        let range = Sphere::new(&Cartessian2D::new([0f64, 0f64]), 1f64);

        let mut a = Cartessian2D::new([-2f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 2f64.sqrt();
        a[1] = 1f64 / 2f64.sqrt();
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        let b = a.clone();
        a[1] = 0f64;
        let movement = Cartessian2D::new([0f64, 3f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        let range = Sphere::new(&CartessianND::new(vec![0f64, 0f64]), 1f64);

        let mut a = CartessianND::new(vec![-2f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 2f64.sqrt();
        a[1] = 1f64 / 2f64.sqrt();
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        let b = a.clone();
        a[1] = 0f64;
        let movement = CartessianND::new(vec![0f64, 3f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );
    }

    #[test]
    fn test_periodic_sphere2d() {
        let range = Sphere::new(&Cartessian2D::new([0f64, 0f64]), 1f64);

        let mut pos = Cartessian2D::new([0f64, 0.9]);
        let mut res = pos.clone();
        assert_abs_diff_eq!(range.find_pair(&pos), &res, epsilon = 1e-3);

        pos[1] = 1.1;
        res[1] = -0.9;
        assert_abs_diff_eq!(range.find_pair(&pos), &res, epsilon = 1e-3);

        pos[0] = 1.0;
        pos[1] = 1.0;

        let sqrt2 = 2f64.sqrt();
        res[0] = 1.0 - sqrt2;
        res[1] = 1.0 - sqrt2;
        assert_abs_diff_eq!(range.find_pair(&pos), &res, epsilon = 1e-3);
    }

    #[test]
    fn test_sphere3d() {
        let range = Sphere::new(&Cartessian3D::new([0f64; 3]), 1f64);

        let mut a = Cartessian3D::new([-2f64, 0f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 2f64.sqrt();
        a[1] = 1f64 / 2f64.sqrt();
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        let b = a.clone();
        a[1] = 0f64;
        let movement = Cartessian3D::new([0f64, 3f64, 0f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        let range = Sphere::new(&CartessianND::new(vec![0f64; 3]), 1f64);

        let mut a = CartessianND::new(vec![-2f64, 0f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 3f64.sqrt();
        a[1] = 1f64 / 3f64.sqrt();
        a[2] = 1f64 / 3f64.sqrt();
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        let b = a.clone();
        a[1] = 0f64;
        let movement = CartessianND::new(vec![0f64, 3f64, 0f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );
    }

    #[test]
    fn test_taxisphere1d() {
        let range = TaxiSphere::new(&Cartessian1D::new([0f64]), 1f64);

        let mut a = Cartessian1D::new([-2f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), Cartessian1D::new([1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), Cartessian1D::new([1f64]));

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), Cartessian1D::new([-1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), Cartessian1D::new([-1f64]));

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 0f64;
        let mut movement = Cartessian1D::new([2f64]);
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            Cartessian1D::new([1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            Cartessian1D::new([1f64])
        );

        movement[0] = -2f64;
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            Cartessian1D::new([-1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            Cartessian1D::new([-1f64])
        );

        let range = TaxiSphere::new(&CartessianND::new(vec![0f64]), 1f64);

        let mut a = CartessianND::new(vec![-2f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), CartessianND::new(vec![1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), CartessianND::new(vec![1f64]));

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), CartessianND::new(vec![-1f64]));
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), CartessianND::new(vec![-1f64]));

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 0f64;
        let mut movement = CartessianND::new(vec![2f64]);
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            CartessianND::new(vec![1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            CartessianND::new(vec![1f64])
        );

        movement[0] = -2f64;
        assert_abs_diff_eq!(
            range.find_intersect(&a, &movement).unwrap(),
            CartessianND::new(vec![-1f64])
        );
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            CartessianND::new(vec![-1f64])
        );
    }

    #[test]
    fn test_taxisphere2d() {
        let range = TaxiSphere::new(&Cartessian2D::new([0f64, 0f64]), 1f64);

        let mut a = Cartessian2D::new([-2f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 3f64;
        a[1] = 2f64 / 3f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(
            range.normal_at(&a).unwrap(),
            Cartessian2D::new([-1f64 / 2f64.sqrt(); 2])
        );
        assert_abs_diff_eq!(
            range.normal_at_unsafe(&a),
            Cartessian2D::new([-1f64 / 2f64.sqrt(); 2])
        );

        let mut b = a.clone();
        a[1] = 0f64;
        let mut movement = Cartessian2D::new([0f64, 3f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        b[1] = -2f64 / 3f64;
        movement[1] = -3f64;
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        movement.coord = [2f64, 0f64];
        b.coord = [1f64, 0f64];
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        let range = TaxiSphere::new(&CartessianND::new(vec![0f64, 0f64]), 1f64);

        let mut a = CartessianND::new(vec![-2f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 3f64;
        a[1] = 2f64 / 3f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(
            range.normal_at(&a).unwrap(),
            CartessianND::new(vec![-1f64 / 2f64.sqrt(); 2])
        );
        assert_abs_diff_eq!(
            range.normal_at_unsafe(&a),
            CartessianND::new(vec![-1f64 / 2f64.sqrt(); 2])
        );

        let mut b = a.clone();
        a[1] = 0f64;
        let mut movement = CartessianND::new(vec![0f64, 3f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        b[1] = -2f64 / 3f64;
        movement[1] = -3f64;
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        movement[0] = 2f64;
        movement[1] = 0f64;
        b[0] = 1f64;
        b[1] = 0f64;
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );
    }

    #[test]
    fn test_taxisphere3d() {
        let range = TaxiSphere::new(&Cartessian3D::new([0f64; 3]), 1f64);

        let mut a = Cartessian3D::new([-2f64, 0f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 3f64;
        a[1] = 1f64 / 3f64;
        a[2] = 1f64 / 3f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(
            range.normal_at(&a).unwrap(),
            Cartessian3D::new([-1f64 / 3f64.sqrt(); 3])
        );
        assert_abs_diff_eq!(
            range.normal_at_unsafe(&a),
            Cartessian3D::new([-1f64 / 3f64.sqrt(); 3])
        );

        a[0] = 1f64 / 3f64;
        a[1] = 2f64 / 3f64;
        a[2] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(
            range.normal_at(&a).unwrap(),
            Cartessian3D::new([-1f64 / 2f64.sqrt(), -1f64 / 2f64.sqrt(), 0f64])
        );
        assert_abs_diff_eq!(
            range.normal_at_unsafe(&a),
            Cartessian3D::new([-1f64 / 2f64.sqrt(), -1f64 / 2f64.sqrt(), 0f64])
        );

        let mut b = a.clone();
        a[1] = 0f64;
        let mut movement = Cartessian3D::new([0f64, 3f64, 0f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        b[1] = -2f64 / 3f64;
        movement[1] = -3f64;
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        movement.coord = [2f64, 0f64, 0f64];
        b.coord = [1f64, 0f64, 0f64];
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        let range = TaxiSphere::new(&CartessianND::new(vec![0f64; 3]), 1f64);

        let mut a = CartessianND::new(vec![-2f64, 0f64, 0f64]);
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = -1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(range.normal_at(&a).unwrap(), -a.clone());
        assert_abs_diff_eq!(range.normal_at_unsafe(&a), -a.clone());

        a[0] = 2f64;
        assert_eq!(range.check_inclusion(&a), false);
        assert_eq!(range.normal_at(&a), None);

        a[0] = 1f64 / 3f64;
        a[1] = 1f64 / 3f64;
        a[2] = 1f64 / 3f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(
            range.normal_at(&a).unwrap(),
            CartessianND::new(vec![-1f64 / 3f64.sqrt(); 3])
        );
        assert_abs_diff_eq!(
            range.normal_at_unsafe(&a),
            CartessianND::new(vec![-1f64 / 3f64.sqrt(); 3])
        );

        a[0] = 1f64 / 3f64;
        a[1] = 2f64 / 3f64;
        a[2] = 0f64;
        assert_eq!(range.check_inclusion(&a), true);
        assert_abs_diff_eq!(
            range.normal_at(&a).unwrap(),
            CartessianND::new(vec![-1f64 / 2f64.sqrt(), -1f64 / 2f64.sqrt(), 0f64])
        );
        assert_abs_diff_eq!(
            range.normal_at_unsafe(&a),
            CartessianND::new(vec![-1f64 / 2f64.sqrt(), -1f64 / 2f64.sqrt(), 0f64])
        );

        let mut b = a.clone();
        a[1] = 0f64;
        let mut movement = CartessianND::new(vec![0f64, 3f64, 0f64]);
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        b[1] = -2f64 / 3f64;
        movement[1] = -3f64;
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );

        movement[0] = 2f64;
        movement[1] = 0f64;
        b[0] = 1f64;
        b[1] = 0f64;
        assert_abs_diff_eq!(range.find_intersect(&a, &movement).unwrap(), b.clone());
        assert_abs_diff_eq!(
            range.find_intersect_unsafe(&a, &movement).unwrap(),
            b.clone()
        );
    }

    #[test]
    fn test_serde_sphere() {
        let center = Cartessian2D::new([0.0, 0.0]);
        let radius = 1.0;
        let sphere: Sphere<Cartessian2D<f64>> = Sphere::new(&center, radius);
        let expected = r#"{"center":{"coord":[0.0,0.0]},"radius":1.0}"#;
        assert_eq!(expected, to_string(&sphere).unwrap());

        let expected = from_str(&expected).unwrap();
        assert_eq!(sphere, expected);
    }
}
