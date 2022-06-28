use crate::force::Bimolecular;
use crate::state::Mass;
use crate::state::State;
use crate::vector::basic::Map;
use crate::vector::product::Distance;
use crate::vector::product::Norm;
use crate::vector::Scalar;
use crate::vector::Vector;
use serde::{Deserialize, Serialize};
use std::ops::AddAssign;
use std::ops::Neg;
use std::ops::Sub;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gravity<T> {
    pub const_g: T,
}

impl<T> PartialEq for Gravity<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.const_g.eq(&other.const_g)
    }
}

impl<T> Gravity<T> {
    pub fn new(const_g: T) -> Self {
        Self { const_g }
    }
}

impl<'a, S, P, T> Bimolecular<'a, S> for Gravity<T>
where
    S: State<Position = P> + Mass<T>,
    P: Vector<Item = T> + Norm<Output = T> + AddAssign<P> + Map<Item = T> + Clone + 'a,
    &'a P: Distance<&'a P, Output = T> + Sub<Output = P> + IntoIterator<Item = &'a T>,
    T: Scalar + Neg<Output = T>,
{
    type Force = P;
    type Potential = T;

    fn potential(&self, state: &'a S, other: &'a S) -> Self::Potential {
        let r = state.pos().distance(other.pos());
        let mass1 = state.mass();
        let mass2 = other.mass();
        return -self.const_g * mass1 * mass2 / r;
    }

    fn force(&self, state: &'a S, other: &'a S) -> Self::Force {
        let r = state.pos().distance(other.pos());
        let c = self.const_g * state.mass() * other.mass() / (r * r * r);
        let mut res = other.pos() - state.pos();
        res.map_inplace(|x| *x = c * *x);

        return res;
    }

    fn force_to(&self, state: &'a S, other: &'a S, force: &'a mut Self::Force) {
        force.clone_from(other.pos());
        force.zip_mut_with(state.pos(), |x, y| *x = *x - *y);
        let r = force.norm_l2();
        let c = self.const_g * state.mass() * other.mass() / (r * r * r);
        force.map_inplace(|x| *x = c * *x);
    }

    fn force_add_to(&self, state: &'a S, other: &'a S, force: &'a mut Self::Force) {
        let mut dist = other.pos() - state.pos();
        let r = dist.norm_l2();
        let c = self.const_g * state.mass() * other.mass() / (r * r * r);
        dist.map_inplace(|x| *x = c * *x);
        *force += dist;
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use serde_json::{from_str, to_string};

    #[test]
    fn test_serde_gravity() {
        let a: Gravity<f64> = Gravity::new(0.5);
        let expected = r#"{"const_g":0.5}"#;
        assert_eq!(expected, to_string(&a).unwrap());

        let expected: Gravity<f64> = from_str(&expected).unwrap();
        assert_eq!(a, expected);
    }
}
