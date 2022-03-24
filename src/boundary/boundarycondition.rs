
#![allow(unused_imports)]
use crate::boundary::IntBoundary;
use crate::state::HasVelocity;
use crate::vector::CartessianND;
use crate::rng::get_uniform;
use crate::vector::product::Dot;
use std::ops::Mul;
use crate::vector::arithmetic::Scalar;
use crate::boundary::Vector;
use crate::boundary::State;
use crate::vector::Cartessian;

use super::NonPeriodic;
use super::{plane::{SimplePlane, Plane, SimplePlanePair, PlanePair, SimpleBox, Cube, Parallelogram}, BoundaryCondition, BCspecies, AfterMove, Boundary, FloatBoundary};

macro_rules! impl_non_periodic_float_boundary_condition {
    ($ty : ident) => {
        impl<B, const N : usize> BoundaryCondition<Cartessian<$ty, N>> for (B, BCspecies<$ty>)
            where B : NonPeriodic + FloatBoundary<Cartessian<$ty, N>>{
            fn check_bc<'a, S>(&self, state : &'a S, movement : &'a mut (Cartessian<$ty, N>, Cartessian<$ty, N>)) -> AfterMove
            where S : State<Movement = (Cartessian<$ty, N>, Cartessian<$ty, N>), Position = Cartessian<$ty, N>> + HasVelocity{
                if let Some(t) = self.0.ratio_to_intersect_unsafe(state.pos(), &movement.0){
                    match &self.1{
                        BCspecies::Absorbing => {
                            return AfterMove::Dead;
                        },
                        BCspecies::PartiallyAbsorbing(p, rng) => {
                            if get_uniform::<$ty>(&mut rng.borrow_mut()) < *p {
                                return AfterMove::Dead;
                            } else {
                                let intersect = state.pos() + &movement.0 * t;
                                let normal = self.0.normal_at_unsafe(&intersect);

                                let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                                movement.0 -= &normal * k;

                                let v = (2.0 as $ty) * (&movement.1 + state.vel()).dot(&normal);
                                movement.1 -= &normal * v;
                                return AfterMove::Survived;
                            }
                        }
                        BCspecies::Reflective => {
                            let intersect = state.pos() + &movement.0 * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                            movement.0 -= &normal * k;

                            let v = (2.0 as $ty) * (&movement.1 + state.vel()).dot(&normal);
                            movement.1 -= &normal * v;
                            return AfterMove::Survived;
                        },
                        BCspecies::InelasticReflection(e) => {
                            let intersect = state.pos() + &movement.0 * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                            movement.0 -= &normal * k;

                            let v = (1 as $ty + e) * (&movement.1 + state.vel()).dot(&normal);
                            movement.1 -= &normal * v;
                            return AfterMove::Survived;
                        },
                        BCspecies::DiffusiveReflection => {
                            unimplemented!();
                        },
                        BCspecies::Periodic => {
                            panic!("Feature is not provided : Periodic boundary condition");
                        },
                    }
                } else {
                    return AfterMove::Survived;
                }
            }

            fn check_bc_overdamped<'a, S>(&self, state : &'a S, movement : &'a mut Cartessian<$ty, N>) -> AfterMove
            where S : State<Movement = Cartessian<$ty, N>, Position = Cartessian<$ty, N>> {
                if let Some(t) = self.0.ratio_to_intersect_unsafe(state.pos(), movement){
                    match &self.1{
                        BCspecies::Absorbing => {
                            return AfterMove::Dead;
                        },
                        BCspecies::PartiallyAbsorbing(p, rng) => {
                            if get_uniform::<$ty>(&mut rng.borrow_mut()) < *p {
                                return AfterMove::Dead;
                            } else {
                                let intersect = state.pos() + movement.clone() * t;
                                let normal = self.0.normal_at_unsafe(&intersect);

                                let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.dot(&normal));
                                *movement -= &normal * k;
                                return AfterMove::Survived;
                            }
                        }
                        BCspecies::Reflective => {
                            let intersect = state.pos() + movement.clone() * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.dot(&normal));
                            *movement -= &normal * k;

                            return AfterMove::Survived;
                        },
                        BCspecies::InelasticReflection(_e) => {
                            panic!("Feature is not provided : Inelastic Reflective boundary condition");
                        },
                        BCspecies::DiffusiveReflection => {
                            unimplemented!();
                        },
                        BCspecies::Periodic => {
                            panic!("Feature is not provided : Periodic boundary condition");
                        },
                    }
                } else {
                    return AfterMove::Survived;
                }
            }
        }

        impl<B> BoundaryCondition<CartessianND<$ty>> for (B, BCspecies<$ty>)
            where B : NonPeriodic + FloatBoundary<CartessianND<$ty>>{
            fn check_bc<'a, S>(&self, state : &'a S, movement : &'a mut (CartessianND<$ty>, CartessianND<$ty>)) -> AfterMove
            where S : State<Movement = (CartessianND<$ty>, CartessianND<$ty>), Position = CartessianND<$ty>> + HasVelocity{
                if let Some(t) = self.0.ratio_to_intersect_unsafe(state.pos(), &movement.0){
                    match &self.1{
                        BCspecies::Absorbing => {
                            return AfterMove::Dead;
                        },
                        BCspecies::PartiallyAbsorbing(p, rng) => {
                            if get_uniform::<$ty>(&mut rng.borrow_mut()) < *p {
                                return AfterMove::Dead;
                            } else {
                                let intersect = state.pos() + &movement.0 * t;
                                let normal = self.0.normal_at_unsafe(&intersect);

                                let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                                movement.0.zip_mut_with(&normal, |x, y| *x = *x - k * *y);

                                let v = (2.0 as $ty) * (&movement.1 + state.vel()).dot(&normal);
                                movement.1.zip_mut_with(&normal, |x, y| *x = *x - v * *y);
                                return AfterMove::Survived;
                            }
                        }
                        BCspecies::Reflective => {
                            let intersect = state.pos() + &movement.0 * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                            movement.0.zip_mut_with(&normal, |x, y| *x = *x - k * *y);

                            let v = (2.0 as $ty) * (&movement.1 + state.vel()).dot(&normal);
                            movement.1.zip_mut_with(&normal, |x, y| *x = *x - v * *y);
                            return AfterMove::Survived;
                        },
                        BCspecies::InelasticReflection(e) => {
                            let intersect = state.pos() + &movement.0 * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                            movement.0.zip_mut_with(&normal, |x, y| *x = *x - k * *y);

                            let v = (1 as $ty + e) * (&movement.1 + state.vel()).dot(&normal);
                            movement.1.zip_mut_with(&normal, |x, y| *x = *x - v * *y);
                            return AfterMove::Survived;
                        },
                        BCspecies::DiffusiveReflection => {
                            unimplemented!();
                        },
                        BCspecies::Periodic => {
                            panic!("Feature is not provided : Periodic boundary condition");
                        },
                    }
                } else {
                    return AfterMove::Survived;
                }
            }

            fn check_bc_overdamped<'a, S>(&self, state : &'a S, movement : &'a mut CartessianND<$ty>) -> AfterMove
            where S : State<Movement = CartessianND<$ty>, Position = CartessianND<$ty>> {
                if let Some(t) = self.0.ratio_to_intersect_unsafe(state.pos(), movement){
                    match &self.1{
                        BCspecies::Absorbing => {
                            return AfterMove::Dead;
                        },
                        BCspecies::PartiallyAbsorbing(p, rng) => {
                            if get_uniform::<$ty>(&mut rng.borrow_mut()) < *p {
                                return AfterMove::Dead;
                            } else {
                                let intersect = state.pos() + movement.clone() * t;
                                let normal = self.0.normal_at_unsafe(&intersect);

                                let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.dot(&normal));
                                movement.zip_mut_with(&normal, |x, y| *x = *x - k * *y);
                                return AfterMove::Survived;
                            }
                        }
                        BCspecies::Reflective => {
                            let intersect = state.pos() + movement.clone() * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k =  (2.0 as $ty) * (1.0 as $ty - t) * (movement.dot(&normal));
                            movement.zip_mut_with(&normal, |x, y| *x = *x - k * *y);

                            return AfterMove::Survived;
                        },
                        BCspecies::InelasticReflection(_e) => {
                            panic!("Feature is not provided : Inelastic Reflective boundary condition");
                        },
                        BCspecies::DiffusiveReflection => {
                            unimplemented!();
                        },
                        BCspecies::Periodic => {
                            panic!("Feature is not provided : Periodic boundary condition");
                        },
                    }
                } else {
                    return AfterMove::Survived;
                }
            }
        }
    };
}

impl_non_periodic_float_boundary_condition!(f32);
impl_non_periodic_float_boundary_condition!(f64);

// impl<B, const N : usize> BoundaryCondition<Cartessian<i32, N>> for (B, BCspecies<i32>)
//             where B : NonPeriodic + IntBoundary<Cartessian<i32, N>>{
//     fn check_bc<'a, S>(&self, state : &'a S, movement : &'a mut (Cartessian<i32, N>, Cartessian<i32, N>)) -> AfterMove
//     where S : State<Movement = (Cartessian<i32, N>, Cartessian<i32, N>), Position = Cartessian<i32, N>> + HasVelocity {

//     }

//     fn check_bc_overdamped<'a, S>(&self, state : &'a S, movement : &'a mut Cartessian<i32, N>) -> AfterMove
//     where S : State<Movement = Cartessian<i32, N>, Position = Cartessian<i32, N>> {

//     }
// }

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;
    use crate::state::Mass;
    use crate::state::HasVelocity;
    use crate::vector::{Cartessian2D, Cartessian3D};
    use crate::prelude::State;
    use crate::boundary::plane::Direction;
    use super::*;

    #[test]
    fn test_bc_simple_plane_float(){
        #[derive(State)]
        struct TestState{
            mass : f64,
            pos : Cartessian2D<f64>,
            vel : Cartessian2D<f64>,
        }

        let state = TestState{
            mass : 1f64,
            pos : Cartessian2D::new([1f64, 1f64]),
            vel : Cartessian2D::new([2f64, -2f64]),
        };

        let plane = SimplePlane::new(1, 0f64, Direction::Positive);

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let mut dv = Cartessian2D::new([0f64, -1f64]);
        let original = (dx.clone(), dv.clone());

        let mut movement = original.clone();
        assert_eq!((plane, BCspecies::Absorbing).check_bc(&state, &mut movement), AfterMove::Dead);

        let mut movement = original.clone();
        assert_eq!((plane, BCspecies::Reflective).check_bc(&state, &mut movement), AfterMove::Survived);
        dx[1] = 0f64; dv[1] = 5f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon=1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon=1e-3);

        let mut movement = original.clone();
        assert_eq!((plane, BCspecies::InelasticReflection(0.5)).check_bc(&state, &mut movement), AfterMove::Survived);
        dx[1] = 0f64; dv[1] = 3.5f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon=1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon=1e-3);

        #[derive(State)]
        struct TestState2{
            mass : f64,
            pos : Cartessian2D<f64>,
        }

        let state = TestState2{
            mass : 1f64,
            pos : Cartessian2D::new([1f64, 1f64]),
        };

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let original = dx.clone();

        let mut movement = original.clone();
        assert_eq!((plane, BCspecies::Absorbing).check_bc_overdamped(&state, &mut movement), AfterMove::Dead);

        let mut movement = original.clone();
        assert_eq!((plane, BCspecies::Reflective).check_bc_overdamped(&state, &mut movement), AfterMove::Survived);
        dx[1] = 0f64;
        assert_abs_diff_eq!(&movement, &dx, epsilon=1e-3);
    }

    #[test]
    fn test_bc_plane_float(){
        #[derive(State)]
        struct TestState{
            mass : f64,
            pos : Cartessian2D<f64>,
            vel : Cartessian2D<f64>,
        }

        let state = TestState{
            mass : 1f64,
            pos : Cartessian2D::new([1f64, 1f64]),
            vel : Cartessian2D::new([2f64, -2f64]),
        };

        let plane = Plane::new(&Cartessian2D::new([0f64, 1f64]), 0f64);

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let mut dv = Cartessian2D::new([0f64, -1f64]);
        let original = (dx.clone(), dv.clone());

        let mut movement = original.clone();
        assert_eq!((plane.clone(), BCspecies::Absorbing).check_bc(&state, &mut movement), AfterMove::Dead);

        let mut movement = original.clone();
        assert_eq!((plane.clone(), BCspecies::Reflective).check_bc(&state, &mut movement), AfterMove::Survived);
        dx[1] = 0f64; dv[1] = 5f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon=1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon=1e-3);

        let mut movement = original.clone();
        assert_eq!((plane.clone(), BCspecies::InelasticReflection(0.5)).check_bc(&state, &mut movement), AfterMove::Survived);
        dx[1] = 0f64; dv[1] = 3.5f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon=1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon=1e-3);

        #[derive(State)]
        struct TestState2{
            mass : f64,
            pos : Cartessian2D<f64>,
        }

        let state = TestState2{
            mass : 1f64,
            pos : Cartessian2D::new([1f64, 1f64]),
        };

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let original = dx.clone();

        let mut movement = original.clone();
        assert_eq!((plane.clone(), BCspecies::Absorbing).check_bc_overdamped(&state, &mut movement), AfterMove::Dead);

        let mut movement = original.clone();
        assert_eq!((plane.clone(), BCspecies::Reflective).check_bc_overdamped(&state, &mut movement), AfterMove::Survived);
        dx[1] = 0f64;
        assert_abs_diff_eq!(&movement, &dx, epsilon=1e-3);
    }

}
