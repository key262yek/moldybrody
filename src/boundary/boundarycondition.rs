#![allow(unused_imports)]
use crate::boundary::sphere::Sphere;
use crate::boundary::IntBoundary;
use crate::boundary::OverdampedBoundaryCondition;
use crate::boundary::Periodic;
use crate::boundary::State;
use crate::boundary::Vector;
use crate::rng::get_uniform;
use crate::state::HasVelocity;
use crate::vector::basic::Map;
use crate::vector::product::Dot;
use crate::vector::Cartessian;
use crate::vector::CartessianND;
use crate::vector::Scalar;
use std::ops::Mul;

use super::NonPeriodic;
use super::{
    plane::{Cube, Parallelogram, Plane, PlanePair, SimpleBox, SimplePlane, SimplePlanePair},
    AfterMove, BCspecies, BoundaryCondition, FloatBoundary,
};

macro_rules! repeat_float_bc {
    ($periodicity : tt, $type : ident<$tok : tt>) => {
        repeat_float_bc!($periodicity, f32, $type<$tok>);
        repeat_float_bc!($periodicity, f64, $type<$tok>);
    };
    ($periodicity : tt, $ty : ident, $type : ident<$tok : tt>) => {
        repeat_float_bc!($periodicity, fixed, $ty, $type<$tok>);
        repeat_float_bc!($periodicity, nd, $ty, $type<$tok>);
    };
    ($periodicity : tt, fixed, $ty : ident, $type : ident<TN>) => {
        impl<S, const N : usize> BoundaryCondition<S, repeat_float_bc!(fixed, $ty)> for ($type<$ty, N>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(fixed, $ty), Movement = (repeat_float_bc!(fixed, $ty), repeat_float_bc!(fixed, $ty))> + HasVelocity,
        {
            repeat_float_bc!($periodicity, fixed, hasvel, $ty);
        }

        impl<S, const N : usize> OverdampedBoundaryCondition<S, repeat_float_bc!(fixed, $ty)> for ($type<$ty, N>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(fixed, $ty), Movement = repeat_float_bc!(fixed, $ty)>
        {
            repeat_float_bc!($periodicity, fixed, overdamped, $ty);
        }
    };
    ($periodicity : tt, nd, $ty : ident, $type : ident<TN>) => {
        impl<S, const N : usize> BoundaryCondition<S, repeat_float_bc!(nd, $ty)> for ($type<$ty, N>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(nd, $ty), Movement = (repeat_float_bc!(nd, $ty), repeat_float_bc!(nd, $ty))> + HasVelocity,
        {
            repeat_float_bc!($periodicity, nd, hasvel, $ty);
        }

        impl<S, const N : usize> OverdampedBoundaryCondition<S, repeat_float_bc!(nd, $ty)> for ($type<$ty, N>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(nd, $ty), Movement = repeat_float_bc!(nd, $ty)>
        {
            repeat_float_bc!($periodicity, nd, overdamped, $ty);
        }
    };
    ($periodicity : tt, fixed, $ty : ident, $type : ident<T>) => {
        impl<S, const N : usize> BoundaryCondition<S, repeat_float_bc!(fixed, $ty)> for ($type<$ty>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(fixed, $ty), Movement = (repeat_float_bc!(fixed, $ty), repeat_float_bc!(fixed, $ty))> + HasVelocity,
        {
            repeat_float_bc!($periodicity, fixed, hasvel, $ty);
        }

        impl<S, const N : usize> OverdampedBoundaryCondition<S, repeat_float_bc!(fixed, $ty)> for ($type<$ty>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(fixed, $ty), Movement = repeat_float_bc!(fixed, $ty)>
        {
            repeat_float_bc!($periodicity, fixed, overdamped, $ty);
        }
    };
    ($periodicity : tt, nd, $ty : ident, $type : ident<T>) => {
        impl<S> BoundaryCondition<S, repeat_float_bc!(nd, $ty)> for ($type<$ty>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(nd, $ty), Movement = (repeat_float_bc!(nd, $ty), repeat_float_bc!(nd, $ty))> + HasVelocity,
        {
            repeat_float_bc!($periodicity, nd, hasvel, $ty);
        }

        impl<S> OverdampedBoundaryCondition<S, repeat_float_bc!(nd, $ty)> for ($type<$ty>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(nd, $ty), Movement = repeat_float_bc!(nd, $ty)>
        {
            repeat_float_bc!($periodicity, nd, overdamped, $ty);
        }
    };
    ($periodicity : tt, fixed, $ty : ident, $type : ident<V>) => {
        impl<S, const N : usize> BoundaryCondition<S, repeat_float_bc!(fixed, $ty)> for ($type<repeat_float_bc!(fixed, $ty)>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(fixed, $ty), Movement = (repeat_float_bc!(fixed, $ty), repeat_float_bc!(fixed, $ty))> + HasVelocity,
        {
            repeat_float_bc!($periodicity, fixed, hasvel, $ty);
        }

        impl<S, const N : usize> OverdampedBoundaryCondition<S, repeat_float_bc!(fixed, $ty)> for ($type<repeat_float_bc!(fixed, $ty)>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(fixed, $ty), Movement = repeat_float_bc!(fixed, $ty)>
        {
            repeat_float_bc!($periodicity, fixed, overdamped, $ty);
        }
    };
    ($periodicity : tt, nd, $ty : ident, $type : ident<V>) => {
        impl<S> BoundaryCondition<S, repeat_float_bc!(nd, $ty)> for ($type<repeat_float_bc!(nd, $ty)>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(nd, $ty), Movement = (repeat_float_bc!(nd, $ty), repeat_float_bc!(nd, $ty))> + HasVelocity,
        {
            repeat_float_bc!($periodicity, nd, hasvel, $ty);
        }

        impl<S> OverdampedBoundaryCondition<S, repeat_float_bc!(nd, $ty)> for ($type<repeat_float_bc!(nd, $ty)>, BCspecies<$ty>)
        where
            S : State<Position = repeat_float_bc!(nd, $ty), Movement = repeat_float_bc!(nd, $ty)>
        {
            repeat_float_bc!($periodicity, nd, overdamped, $ty);
        }
    };
    ($periodicity : tt, $vector : tt, hasvel, $ty : ident) => {
        fn check_bc<'a>(
            &self,
            state: &'a S,
            movement: &'a mut (repeat_float_bc!($vector, $ty), repeat_float_bc!($vector, $ty)),
        ) -> AfterMove
        {
            if let Some(t) = self.0.ratio_to_intersect_unsafe(state.pos(), &movement.0) {
                match &self.1 {
                    BCspecies::Absorbing => {
                        return AfterMove::Dead;
                    }
                    BCspecies::PartiallyAbsorbing(p, rng) => {
                        if get_uniform::<$ty>(&mut rng.write().unwrap()) < *p {
                            return AfterMove::Dead;
                        } else {
                            let intersect = state.pos() + &movement.0 * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k = (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                            movement.0 -= &normal * k;

                            let v = (2.0 as $ty) * (&movement.1 * t + state.vel()).dot(&normal);
                            movement.1 -= &normal * v;
                            return AfterMove::Survived;
                        }
                    }
                    BCspecies::Reflective => {
                        let intersect = state.pos() + &movement.0 * t;
                        let normal = self.0.normal_at_unsafe(&intersect);

                        let k = (2.0 as $ty) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                        movement.0 -= &normal * k;

                        let v = (2.0 as $ty) * (&movement.1 * t + state.vel()).dot(&normal);
                        movement.1 -= &normal * v;
                        return AfterMove::Survived;
                    }
                    BCspecies::InelasticReflection(e) => {
                        let intersect = state.pos() + &movement.0 * t;
                        let normal = self.0.normal_at_unsafe(&intersect);

                        let k = (1 as $ty + e) * (1.0 as $ty - t) * (movement.0.dot(&normal));
                        movement.0 -= &normal * k;

                        let v = (1 as $ty + e) * (&movement.1 * t + state.vel()).dot(&normal);
                        movement.1 -= &normal * v;
                        return AfterMove::Survived;
                    }
                    BCspecies::DiffusiveReflection => {
                        unimplemented!();
                    }
                    BCspecies::Periodic => {
                        repeat_float_bc!($periodicity, hasvel, self, state, movement);
                    }
                }
            } else {
                return AfterMove::Survived;
            }
        }
    };
    ($periodicity : tt, $vector : tt, overdamped, $ty : ident) => {
        fn check_bc_overdamped<'a>(
            &self,
            state: &'a S,
            movement: &'a mut repeat_float_bc!($vector, $ty),
        ) -> AfterMove
        {
            if let Some(t) = self.0.ratio_to_intersect_unsafe(state.pos(), movement) {
                match &self.1 {
                    BCspecies::Absorbing => {
                        return AfterMove::Dead;
                    }
                    BCspecies::PartiallyAbsorbing(p, rng) => {
                        if get_uniform::<$ty>(&mut rng.write().unwrap()) < *p {
                            return AfterMove::Dead;
                        } else {
                            let intersect = state.pos() + movement.clone() * t;
                            let normal = self.0.normal_at_unsafe(&intersect);

                            let k = (2.0 as $ty) * (1.0 as $ty - t) * (movement.dot(&normal));
                            *movement -= &normal * k;
                            return AfterMove::Survived;
                        }
                    }
                    BCspecies::Reflective => {
                        let intersect = state.pos() + movement.clone() * t;
                        let normal = self.0.normal_at_unsafe(&intersect);

                        let k = (2.0 as $ty) * (1.0 as $ty - t) * (movement.dot(&normal));
                        *movement -= &normal * k;

                        return AfterMove::Survived;
                    }
                    BCspecies::InelasticReflection(_e) => {
                        panic!("Feature is not provided : Inelastic Reflective boundary condition");
                    }
                    BCspecies::DiffusiveReflection => {
                        unimplemented!();
                    }
                    BCspecies::Periodic => {
                        repeat_float_bc!($periodicity, overdamped, self, state, movement);
                    }
                }
            } else {
                return AfterMove::Survived;
            }
        }
    };
    (fixed, $ty : ident) => {
        Cartessian<$ty, N>
    };
    (nd, $ty : ident) => {
        CartessianND<$ty>
    };
    (periodic, hasvel, $self : ident, $state : ident, $movement : ident) => {
        let mut arrive = $state.pos() + &$movement.0;
        $self.0.find_pair_mut(&mut arrive);
        arrive -= $state.pos();
        $movement.0.clone_from(&arrive);
        return AfterMove::Survived
    };
    (periodic, overdamped, $self : ident, $state : ident, $movement : ident) => {
        let mut arrive = $state.pos() + &*$movement;
        $self.0.find_pair_mut(&mut arrive);
        arrive -= $state.pos();
        $movement.clone_from(&arrive);
        return AfterMove::Survived
    };
    (non_periodic, hasvel, $self : ident, $state : ident, $movement : ident) => {
        panic!("Feature is not provided : Periodic boundary condition")
    };
    (non_periodic, overdamped, $self : ident, $state : ident, $movement : ident) => {
        panic!("Feature is not provided : Periodic boundary condition")
    };
}

repeat_float_bc!(non_periodic, SimplePlane<T>);
repeat_float_bc!(non_periodic, Plane<V>);
repeat_float_bc!(periodic, SimplePlanePair<T>);
repeat_float_bc!(periodic, PlanePair<V>);
repeat_float_bc!(periodic, SimpleBox<TN>);
repeat_float_bc!(periodic, Cube<V>);
repeat_float_bc!(periodic, Sphere<V>);
repeat_float_bc!(periodic, fixed, f32, Parallelogram<TN>);
repeat_float_bc!(periodic, fixed, f64, Parallelogram<TN>);

#[cfg(test)]
mod test {
    use super::*;
    use crate::boundary::plane::Direction;
    use crate::prelude::State;
    use crate::state::HasVelocity;
    use crate::state::Mass;
    use crate::vector::{Cartessian2D, Cartessian3D};
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_bc_simple_plane_float() {
        #[derive(State)]
        struct TestState {
            mass: f64,
            pos: Cartessian2D<f64>,
            vel: Cartessian2D<f64>,
        }

        let state = TestState {
            mass: 1f64,
            pos: Cartessian2D::new([1f64, 1f64]),
            vel: Cartessian2D::new([2f64, -2f64]),
        };

        let plane = SimplePlane::new(1, 0f64, Direction::Positive);

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let mut dv = Cartessian2D::new([0f64, -1f64]);
        let original = (dx.clone(), dv.clone());

        let mut movement = original.clone();
        assert_eq!(
            (plane, BCspecies::Absorbing).check_bc(&state, &mut movement),
            AfterMove::Dead
        );

        let mut movement = original.clone();
        assert_eq!(
            (plane, BCspecies::Reflective).check_bc(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = 0f64;
        dv[1] = 4f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon = 1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon = 1e-3);

        let mut movement = original.clone();
        assert_eq!(
            (plane, BCspecies::InelasticReflection(0.5)).check_bc(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = -0.5f64;
        dv[1] = 2.75f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon = 1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon = 1e-3);

        #[derive(State)]
        struct TestState2 {
            mass: f64,
            pos: Cartessian2D<f64>,
        }

        let state = TestState2 {
            mass: 1f64,
            pos: Cartessian2D::new([1f64, 1f64]),
        };

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let original = dx.clone();

        let mut movement = original.clone();
        assert_eq!(
            (plane, BCspecies::Absorbing).check_bc_overdamped(&state, &mut movement),
            AfterMove::Dead
        );

        let mut movement = original.clone();
        assert_eq!(
            (plane, BCspecies::Reflective).check_bc_overdamped(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = 0f64;
        assert_abs_diff_eq!(&movement, &dx, epsilon = 1e-3);
    }

    #[test]
    fn test_bc_plane_float() {
        #[derive(State)]
        struct TestState {
            mass: f64,
            pos: Cartessian2D<f64>,
            vel: Cartessian2D<f64>,
        }

        let state = TestState {
            mass: 1f64,
            pos: Cartessian2D::new([1f64, 1f64]),
            vel: Cartessian2D::new([2f64, -2f64]),
        };

        let plane = Plane::new(Cartessian2D::new([0f64, 1f64]), 0f64);

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let mut dv = Cartessian2D::new([0f64, -1f64]);
        let original = (dx.clone(), dv.clone());

        let mut movement = original.clone();
        assert_eq!(
            (plane.clone(), BCspecies::Absorbing).check_bc(&state, &mut movement),
            AfterMove::Dead
        );

        let mut movement = original.clone();
        assert_eq!(
            (plane.clone(), BCspecies::Reflective).check_bc(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = 0f64;
        dv[1] = 4f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon = 1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon = 1e-3);

        let mut movement = original.clone();
        assert_eq!(
            (plane.clone(), BCspecies::InelasticReflection(0.5)).check_bc(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = -0.5f64;
        dv[1] = 2.75f64;
        assert_abs_diff_eq!(&movement.0, &dx, epsilon = 1e-3);
        assert_abs_diff_eq!(&movement.1, &dv, epsilon = 1e-3);

        #[derive(State)]
        struct TestState2 {
            mass: f64,
            pos: Cartessian2D<f64>,
        }

        let state = TestState2 {
            mass: 1f64,
            pos: Cartessian2D::new([1f64, 1f64]),
        };

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let original = dx.clone();

        let mut movement = original.clone();
        assert_eq!(
            (plane.clone(), BCspecies::Absorbing).check_bc_overdamped(&state, &mut movement),
            AfterMove::Dead
        );

        let mut movement = original.clone();
        assert_eq!(
            (plane.clone(), BCspecies::Reflective).check_bc_overdamped(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = 0f64;
        assert_abs_diff_eq!(&movement, &dx, epsilon = 1e-3);
    }

    #[test]
    fn test_bc_periodic() {
        #[derive(State)]
        struct TestState {
            mass: f64,
            pos: Cartessian2D<f64>,
        }

        let state = TestState {
            mass: 1f64,
            pos: Cartessian2D::new([1f64, 1f64]),
        };

        let planepair = SimplePlanePair::new(1, [0f64, 2f64]).unwrap();

        let mut dx = Cartessian2D::new([1f64, -2f64]);
        let original = dx.clone();

        let mut movement = original.clone();
        assert_eq!(
            (planepair.clone(), BCspecies::Absorbing).check_bc_overdamped(&state, &mut movement),
            AfterMove::Dead
        );

        let mut movement = original.clone();
        assert_eq!(
            (planepair.clone(), BCspecies::Reflective).check_bc_overdamped(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = 0f64;
        assert_abs_diff_eq!(&movement, &dx, epsilon = 1e-3);

        let mut movement = original.clone();
        assert_eq!(
            (planepair.clone(), BCspecies::Periodic).check_bc_overdamped(&state, &mut movement),
            AfterMove::Survived
        );
        dx[1] = 0f64;
        assert_abs_diff_eq!(&movement, &dx, epsilon = 1e-3);
    }
}
