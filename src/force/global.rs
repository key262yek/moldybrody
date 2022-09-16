use super::Global;
use crate::argument::CommandBuilder;
use crate::force::Vector;
use crate::state::Charge;
use crate::state::HasVelocity;
use crate::state::Mass;
use crate::state::State;
use crate::vector::product::Cross;
use crate::vector::product::Dot;
use crate::vector::Dim;
use crate::vector::Scalar;
use clap::Arg;
use clap::ArgMatches;
use core::str::FromStr;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;
use std::ops::Neg;
use crate::vector::basic::Map;


#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConstGravity<V: Vector> {
    pub acc: V,
}

impl<V> PartialEq for ConstGravity<V>
where
    V: Vector + PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.acc.eq(&other.acc)
    }
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
    V: Vector<Item = T> + Dot<&'a V, Output = T> + Mul<T, Output = V> + Map<Item = T> + Index<usize, Output = T> + Clone + 'a + Debug,
    &'a mut V: IntoIterator<Item = &'a mut T>,
    T: Scalar + Neg<Output = T>,
{
    type Force = V;
    type Potential = T;

    fn potential(&self, state: &'a S) -> Self::Potential {
        -state.mass() * self.acc.dot(state.pos())
    }

    fn force(&self, state: &'a S) -> Self::Force {
        let mut res = self.acc.clone();
        let m = state.mass();
        res.map_inplace(|x| *x *= m);
        res
    }

    fn force_to(&self, state: &'a S, force: &'a mut Self::Force) {
        let m = state.mass();
        for (i, f) in force.into_iter().enumerate(){
            *f = m * self.acc[i];
        }
    }

    fn force_add_to(&self, state: &'a S, force: &'a mut Self::Force) {
        let m = state.mass();
        for (i, f) in force.into_iter().enumerate(){
            *f += m * self.acc[i];
        }
    }
}

impl<'h, V> CommandBuilder<'h, 1> for ConstGravity<V>
where
    V: Vector,
{
    fn args() -> [Arg<'h>; 1] {
        [Arg::new("acc")
            .short('g')
            .long("acc")
            .value_name("ACC")
            .takes_value(true)
            .help("Gravitational Acceleration")]
    }
}

impl<V> From<&ArgMatches> for ConstGravity<V>
where
    V: Vector + FromStr + Clone + Send + Sync + 'static,
    <V as FromStr>::Err: Debug,
{
    fn from(m: &ArgMatches) -> Self {
        let acc: V = m.get_one::<String>("acc").unwrap().parse().unwrap();
        ConstGravity::<V>::new(acc)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Lorentz<V: Vector + Dim<3>> {
    pub mag_field: V,
}

impl<V: Vector + Dim<3>> Lorentz<V> {
    #[allow(dead_code)]
    pub fn new(mag_field: V) -> Self {
        Self { mag_field }
    }
}

impl<V> PartialEq for Lorentz<V>
where
    V: Vector + Dim<3> + PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.mag_field.eq(&other.mag_field)
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
        + 'a + Debug,
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

    fn force_to(& self, state: &'a S, force: &mut Self::Force) {
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

    fn force_add_to(& self, state: &'a S, force: &mut Self::Force) {
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

impl<'h, V> CommandBuilder<'h, 1> for Lorentz<V>
where
    V: Vector + Dim<3>,
{
    fn args() -> [Arg<'h>; 1] {
        [Arg::new("mag_field")
            .short('B')
            .long("mag")
            .value_name("BFIELD")
            .takes_value(true)
            .help("Uniform magnetic field")]
    }
}

impl<V> From<&ArgMatches> for Lorentz<V>
where
    V: Vector + Dim<3> + FromStr + Clone + Send + Sync + 'static,
    <V as FromStr>::Err: Debug,
{
    fn from(m: &ArgMatches) -> Self {
        let mag_field: V = m.get_one::<String>("mag_field").unwrap().parse().unwrap();
        Lorentz::<V>::new(mag_field)
    }
}


struct GeneralGlobal<S: State, F: Vector, P> {
    potential: fn(&S) -> P,
    force: fn(&S) -> F,
}

impl<S, F, P> Debug for GeneralGlobal<S, F, P>
    where S : State,
    F : Vector,{
    fn fmt(&self, f : &mut std::fmt::Formatter) -> std::fmt::Result{
        write!(f, "General Global interaction")
    }
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

    fn potential(&self, state: &'a S) -> Self::Potential {
        (self.potential)(state)
    }

    fn force(&self, state: &'a S) -> Self::Force {
        (self.force)(state)
    }

    fn force_to(&self, state: &'a S, force: &mut Self::Force) {
        *force = (self.force)(state).clone();
    }

    fn force_add_to(&self, state: &'a S, force: &mut Self::Force) {
        *force += (self.force)(state);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::basic::Map;
    use crate::vector::Cartessian;
    use crate::vector::{Cartessian2D, Cartessian3D};
    use approx::assert_abs_diff_eq;
    use clap::builder::Command;
    use moldybrody_proc::State;
    use serde_json::from_str;
    use serde_json::to_string;

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
    fn test_serde_const_gravity() {
        let acc = Cartessian::new([0f64, 10f64]);
        let a: ConstGravity<Cartessian<f64, 2>> = ConstGravity::new(acc);
        let expected = r#"{"acc":{"coord":[0.0,10.0]}}"#;
        assert_eq!(expected, to_string(&a).unwrap());

        let expected: ConstGravity<Cartessian<f64, 2>> = from_str(&expected).unwrap();
        assert_eq!(a, expected);
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
        #[derive(State, Debug)]
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

    #[test]
    fn test_serde_lorentz() {
        let mag_field = Cartessian::new([0f64, 0f64, 10f64]);
        let a: Lorentz<Cartessian<f64, 3>> = Lorentz::new(mag_field);
        let expected = r#"{"mag_field":{"coord":[0.0,0.0,10.0]}}"#;
        assert_eq!(expected, to_string(&a).unwrap());

        let expected: Lorentz<Cartessian<f64, 3>> = from_str(&expected).unwrap();
        assert_eq!(a, expected);
    }

    #[test]
    fn test_clap_global() {
        let arg = Command::new("test")
            .args(ConstGravity::<Cartessian2D<f64>>::args())
            .get_matches_from(vec!["test", "--acc", "0,-10"]);
        let gravity = ConstGravity::<Cartessian2D<f64>>::from(&arg);
        assert_abs_diff_eq!(gravity.acc, Cartessian2D::<f64>::new([0.0, -10.0]));

        let arg = Command::new("test2")
            .args(Lorentz::<Cartessian3D<f64>>::args())
            .get_matches_from(vec!["test2", "--mag", "0,0,-10"]);
        let lorentz = Lorentz::<Cartessian3D<f64>>::from(&arg);
        assert_abs_diff_eq!(
            lorentz.mag_field,
            Cartessian3D::<f64>::new([0.0, 0.0, -10.0])
        );
    }
}
