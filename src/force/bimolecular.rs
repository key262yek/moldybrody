
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::AddAssign;
use std::ops::SubAssign;
use std::ops::Sub;
use std::ops::Neg;
use crate::vector::Vector;
use crate::state::State;
use crate::force::Bimolecular;
use crate::state::Mass;
use crate::vector::product::Norm;
use crate::vector::product::Distance;
use crate::vector::arithmetic::Scalar;

pub struct Gravity<T>{
	const_g : T,
}

impl<T> Gravity<T>{
	pub fn new(const_g : T) -> Self{
		Self{
			const_g,
		}
	}
}

impl<'a, S, M, P, T> Bimolecular<'a, S> for Gravity<T>
	where S : State<Movement = M, Position = P> + Mass<T>,
		  M : Vector<Item = T>,
		  P : Vector<Item = T> + Norm<Output = T> + MulAssign<T> + AddAssign<&'a P> + SubAssign<&'a P> + 'a,
		  &'a P : Distance<&'a P, Output = T> + Sub<Output = P>,
		  T : Scalar + Neg<Output = T> + Mul<P, Output = P>{
	type Force = P;
	type Potential = T; 

	fn potential(&self, state : &'a S, other : &'a S) -> Self::Potential{
		let r = state.pos().distance(other.pos());
		let mass1 = state.mass();
		let mass2 = other.mass();
		return -self.const_g * mass1 * mass2 / r;
	}

    fn force(&self, state : &'a S, other : &'a S) -> Self::Force{
    	let r = state.pos().distance(other.pos());
		let mass1 = state.mass();
		let mass2 = other.mass();
		return self.const_g * mass1 * mass2 / (r * r * r) * (other.pos() - state.pos());	
    }

    
    fn force_to(&self, state : &'a S, other : &'a S, force : &'a mut Self::Force){
    	*force *= T::zero();
    	*force -= state.pos();
    	*force += other.pos();
    	let r = force.norm_l2();
    	*force *= self.const_g * state.mass() * other.mass() / (r * r * r);
    }
}