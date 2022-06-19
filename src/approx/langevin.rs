use crate::approx::ApproxOverdampedLangevin;
use crate::approx::State;
use crate::approx::Vector;
use crate::state::DiffusionFloat;
use std::ops::Add;
use std::ops::Mul;

macro_rules! impl_langevin_float {
    ($ty : ident) => {
        impl<'a, S, P> ApproxOverdampedLangevin<'a, P, $ty> for S
        where
            S: State<Movement = P, Position = P> + DiffusionFloat<$ty>,
            P: Vector<Item = $ty> + Add<Output = P> + 'a,
            &'a P: Mul<$ty, Output = P>,
        {
            /// Pure diffusion
            fn euler_pure_diffusion(&'a self, random_force: &'a P, dt: $ty) -> P {
                random_force * (self.diff_const() * dt.sqrt())
            }

            fn euler_with_drift(&'a self, force: &'a P, random_force: &'a P, dt: $ty) -> P {
                force * (self.diff_const() * dt) + random_force * (self.diff_const() * dt.sqrt())
            }
        }
    };
}

impl_langevin_float!(f32);
impl_langevin_float!(f64);
