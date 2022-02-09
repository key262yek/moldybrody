
use std::ops::AddAssign;
use std::ops::MulAssign;
use std::ops::Add;
use crate::prelude::Mul;
use crate::vector::arithmetic::Scalar;
use crate::state::Mass;
use crate::approx::Vector;
use crate::approx::State;
use crate::state::HasVelocity;
use super::ApproxNewton;



impl<'a, S, P, T> ApproxNewton<'a> for S
    where S : State<Movement = (P, P), Position = P> + HasVelocity + Mass<T>,
          P : Vector<Item = T> + Add<Output = P> + Clone + AddAssign<&'a P> + AddAssign<P>+ MulAssign<T>  + 'a,
          &'a P : Mul<T, Output = P> ,
          T : Scalar{

    fn approx(&'a self, force : &'a P, dt : T) -> (P, P){
        let dv = force * (dt / self.mass());
        let dx = self.vel() * dt + force * (dt * dt / (self.mass() + self.mass()));
        return (dx, dv);
    }
}
