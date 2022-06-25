use super::Global;
use crate::force::Vector;
use crate::state::Charge;
use crate::state::HasVelocity;
use crate::state::Mass;
use crate::state::State;
use crate::vector::product::Cross;
use crate::vector::product::Dot;
use crate::vector::Dim;
use crate::vector::Scalar;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;
use std::ops::Neg;

#[derive(Clone)]
pub struct ConstGravity<V: Vector> {
    pub acc: V,
}

impl<V: Vector> ConstGravity<V> {
    #[allow(dead_code)]
    pub fn new(acc: V) -> Self {
        Self { acc }
    }
}

impl<'a, S, V, T> Global<'a, S> for ConstGravity<V>
where
    S: State<Position = V> + Mass<T>,
    V: Vector<Item = T> + Dot<&'a V, Output = T> + 'a,
    &'a V: Mul<T, Output = V> + IntoIterator<Item = &'a T>,
    &'a mut V: IntoIterator<Item = &'a mut T>,
    T: Scalar + Neg<Output = T>,
{
    type Force = V;
    type Potential = T;

    fn potential(&'a self, state: &'a S) -> Self::Potential {
        -state.mass() * self.acc.dot(state.pos())
    }

    fn force(&'a self, state: &'a S) -> Self::Force {
        &self.acc * state.mass()
    }

    fn force_to(&'a self, state: &'a S, force: &'a mut Self::Force) {
        for (f, s) in force.into_iter().zip(self.acc.into_iter()) {
            *f = state.mass() * *s;
        }
    }

    fn force_add_to(&'a self, state: &'a S, force: &'a mut Self::Force) {
        for (f, s) in force.into_iter().zip(self.acc.into_iter()) {
            *f += state.mass() * *s;
        }
    }
}

pub struct Lorentz<V: Vector + Dim<3>> {
    pub mag_field: V,
}

impl<V: Vector + Dim<3>> Lorentz<V> {
    #[allow(dead_code)]
    pub fn new(mag_field: V) -> Self {
        Self { mag_field }
    }
}

impl<'a, S, V, T> Global<'a, S> for Lorentz<V>
where
    S: State<Position = V> + Charge<T> + HasVelocity,
    V: Vector<Item = T>
        + Dim<3>
        + Cross<&'a V, Output = V>
        + Index<usize, Output = T>
        + IndexMut<usize>
        + Div<T, Output = V>
        + 'a,
    T: Scalar + Mul<V, Output = V> + Neg<Output = T>,
{
    type Force = V;
    type Potential = V;

    fn potential(&self, state: &'a S) -> Self::Potential {
        // calculate vector potential
        let two = T::one() + T::one();
        return self.mag_field.cross(&state.pos()) / two;
    }

    fn force(&self, state: &'a S) -> Self::Force {
        -state.charge() * self.mag_field.cross(&state.vel())
    }

    fn force_to(&'a self, state: &'a S, force: &'a mut Self::Force) {
        if self.mag_field.dim() != state.vel().dim() {
            panic!("Cross product between different length of vector is not available.");
        } else if self.mag_field.dim() != 3 {
            panic!("Lorentz force is available only for 3D");
        }

        let q = state.charge();
        let vel = state.vel();
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            force[k] = q * (vel[i] * self.mag_field[j] - vel[j] * self.mag_field[i]);
        }
    }

    fn force_add_to(&'a self, state: &'a S, force: &'a mut Self::Force) {
        if self.mag_field.dim() != state.vel().dim() {
            panic!("Cross product between different length of vector is not available.");
        } else if self.mag_field.dim() != 3 {
            panic!("Lorentz force is available only for 3D");
        }

        let q = state.charge();
        let vel = state.vel();
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            force[k] += q * (vel[i] * self.mag_field[j] - vel[j] * self.mag_field[i]);
        }
    }
}

struct GeneralGlobal<S: State, F: Vector, P> {
    potential: fn(&S) -> P,
    force: fn(&S) -> F,
}

impl<S, F, P> GeneralGlobal<S, F, P>
where
    S: State,
    F: Vector,
{
    #[allow(dead_code)]
    pub fn new(potential: fn(&S) -> P, force: fn(&S) -> F) -> Self {
        GeneralGlobal { potential, force }
    }
}

impl<'a, S, F, P> Global<'a, S> for GeneralGlobal<S, F, P>
where
    S: State,
    F: Vector + AddAssign + Clone,
{
    type Force = F;
    type Potential = P;

    fn potential(&'a self, state: &'a S) -> Self::Potential {
        (self.potential)(state)
    }

    fn force(&'a self, state: &'a S) -> Self::Force {
        (self.force)(state)
    }

    fn force_to(&'a self, state: &'a S, force: &'a mut Self::Force) {
        *force = (self.force)(state).clone();
    }

    fn force_add_to(&'a self, state: &'a S, force: &'a mut Self::Force) {
        *force += (self.force)(state);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::basic::Map;
    use crate::vector::{Cartessian2D, Cartessian3D};
    use approx::assert_abs_diff_eq;
    use moldybrody_proc::State;

    #[test]
    fn test_constgravity() {
        struct TestState {
            mass: f64,
            pos: Cartessian2D<f64>,
        }

        impl Mass<f64> for TestState {
            fn mass(&self) -> f64 {
                self.mass
            }
        }

        impl State for TestState {
            type Movement = Cartessian2D<f64>;
            type Position = Cartessian2D<f64>;

            fn pos(&self) -> &Self::Position {
                &self.pos
            }

            fn pos_mut(&mut self) -> &mut Self::Position {
                &mut self.pos
            }

            fn disp<'a>(&self, movement: &'a Self::Movement) -> &'a Self::Position {
                movement
            }

            fn renew_state(&mut self, movement: &Self::Movement) {
                self.pos.zip_mut_with(movement, |p, m| *p = *p + *m);
            }
        }

        let gravity = ConstGravity::new(Cartessian2D::new([0f64, -10f64]));
        let state = TestState {
            mass: 10f64,
            pos: Cartessian2D::new([5f64, 10f64]),
        };
        let mut force = Cartessian2D::<f64>::default();

        assert_eq!(gravity.potential(&state), 1000f64);
        assert_eq!(gravity.force(&state), Cartessian2D::new([0f64, -100f64]));
        gravity.force_to(&state, &mut force);
        assert_eq!(force, Cartessian2D::new([0f64, -100f64]));
    }

    #[test]
    fn test_lorentz() {
        struct TestState {
            #[allow(dead_code)]
            mass: f64,
            charge: f64,
            pos: Cartessian3D<f64>,
            vel: Cartessian3D<f64>,
        }

        impl Charge<f64> for TestState {
            fn charge(&self) -> f64 {
                self.charge
            }
        }

        impl State for TestState {
            type Movement = (Cartessian3D<f64>, Cartessian3D<f64>);
            type Position = Cartessian3D<f64>;

            fn pos(&self) -> &Self::Position {
                &self.pos
            }

            fn pos_mut(&mut self) -> &mut Self::Position {
                &mut self.pos
            }

            fn disp<'a>(&self, movement: &'a Self::Movement) -> &'a Self::Position {
                &movement.0
            }

            fn renew_state(&mut self, movement: &Self::Movement) {
                self.pos.zip_mut_with(&movement.0, |p, m| *p = *p + *m);
                self.vel.zip_mut_with(&movement.1, |v, m| *v = *v + *m);
            }
        }

        impl HasVelocity for TestState {
            fn vel(&self) -> &<Self as State>::Position {
                &self.vel
            }

            fn vel_mut(&mut self) -> &mut <Self as State>::Position {
                &mut self.vel
            }
        }

        let lorentz = Lorentz::new(Cartessian3D::new([0f64, 0f64, 10f64]));
        let state = TestState {
            mass: 1e-1f64,
            charge: 1e-3f64,
            pos: Cartessian3D::new([5f64, 0f64, 0f64]),
            vel: Cartessian3D::new([0f64, -0.5f64, 0f64]),
        };
        let mut force = Cartessian3D::<f64>::default();

        assert_eq!(
            lorentz.potential(&state),
            Cartessian3D::new([0f64, 25f64, 0f64])
        );
        assert_eq!(
            lorentz.force(&state),
            Cartessian3D::new([-5e-3f64, 0f64, 0f64])
        );
        lorentz.force_to(&state, &mut force);
        assert_eq!(force, Cartessian3D::new([-5e-3f64, 0f64, 0f64]));
    }

    #[test]
    fn general_global_interaction() {
        #[derive(State)]
        struct TestState {
            mass: f64,
            pos: Cartessian2D<f64>,
        }

        fn pot(state: &TestState) -> f64 {
            state.mass() * state.pos()[1] * 10.0
        }

        fn force(state: &TestState) -> Cartessian2D<f64> {
            let acc: Cartessian2D<f64> = Cartessian2D::new([0f64, -10f64]);
            state.mass() * acc
        }

        let interaction: GeneralGlobal<TestState, Cartessian2D<f64>, f64> =
            GeneralGlobal::new(pot, force);
        let state = TestState {
            mass: 1.0,
            pos: Cartessian2D::new([0f64, 0f64]),
        };

        assert_abs_diff_eq!(interaction.potential(&state), 0f64);
        assert_abs_diff_eq!(interaction.force(&state), Cartessian2D::new([0f64, -10f64]));
    }
}
