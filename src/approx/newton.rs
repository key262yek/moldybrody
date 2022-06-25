use crate::approx::ApproxNewton;
use crate::prelude::Mul;
use crate::vector::basic::Zeros;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::MulAssign;

use crate::approx::State;
use crate::approx::Vector;
use crate::state::HasVelocity;
use crate::state::Mass;

macro_rules! impl_newton_float {
    ($ty : ident) => {
        impl<'a, S, P> ApproxNewton<'a, P, $ty> for S
        where
            S: State<Movement = (P, P), Position = P> + HasVelocity + Mass<$ty>,
            P: Vector<Item = $ty>
                + Add<Output = P>
                + Clone
                + AddAssign<&'a P>
                + AddAssign<P>
                + MulAssign<$ty>
                + Zeros
                + Mul<$ty, Output = P>
                + 'a,
            &'a P: Mul<$ty, Output = P> + Add<P, Output = P>,
        {
            fn euler(&'a self, force: &'a P, dt: $ty) -> (P, P) {
                let dv = force * (dt / self.mass());
                let dx = (self.vel() + dv.clone() * 0.5 as $ty) * dt;
                return (dx, dv);
            }

            fn const_speed(&'a self, dt: $ty) -> (P, P) {
                let dx = self.vel() * dt;
                let dv = P::zero_with_length(self.vel().dim());
                return (dx, dv);
            }
        }
    };
}

impl_newton_float!(f32);
impl_newton_float!(f64);
