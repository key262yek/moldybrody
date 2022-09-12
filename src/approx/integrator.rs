use crate::approx::tableau::ApproxMethod;
use crate::approx::tableau::ButcherTableau;
use crate::approx::TimeIterator;
use crate::boundary::System;
use crate::force::{Bimolecular, Global, RandomForce};
use crate::state::State;
use crate::vector::Scalar;
use crate::vector::Vector;
use num_traits::Zero;
use std::ops::Mul;

pub struct RKIntegratorBuilder<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    mode: Option<&'a str>,
    method: Option<Box<dyn ApproxMethod<<V as Vector>::Item>>>,
    time_iter: Option<T>,
    states: Option<Vec<S>>,
    global_force: Option<Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    bimolecular_force:
        Option<Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    random_force: Option<Box<dyn RandomForce<'a, S, Force = V>>>,
    boundary_cond: Option<System<S, V>>,
}
impl<'a, V, S, T> RKIntegratorBuilder<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    pub fn new() -> Self {
        Self {
            mode: None,
            method: None,
            time_iter: None,
            states: None,
            global_force: None,
            bimolecular_force: None,
            random_force: None,
            boundary_cond: None,
        }
    }

    pub fn mode(mut self, mode: &'a str) -> Self {
        self.mode = Some(mode);
        self
    }

    pub fn method(mut self, m: Box<dyn ApproxMethod<<V as Vector>::Item>>) -> Self {
        self.method = Some(m);
        self
    }

    pub fn time_iter(mut self, t: T) -> Self {
        self.time_iter = Some(t);
        self
    }

    pub fn states(mut self, states: Vec<S>) -> Self {
        self.states = Some(states);
        self
    }

    pub fn global_force(
        mut self,
        force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    ) -> Self {
        self.global_force = Some(force);
        self
    }

    pub fn bimolecular_force(
        mut self,
        force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    ) -> Self {
        self.bimolecular_force = Some(force);
        self
    }

    pub fn random_force(mut self, force: Box<dyn RandomForce<'a, S, Force = V>>) -> Self {
        self.random_force = Some(force);
        self
    }

    pub fn boundary_cond(mut self, sys: System<S, V>) -> Self {
        self.boundary_cond = Some(sys);
        self
    }

    pub fn build(self) -> RKIntegrator<'a, V, S, T> {
        // match self.mode {
        // 	Some("Newton") => {
        // 		match (self.global_force, self.bimolecular_force) {
        // 			(Some(g), Some(b)) =>
        // 			(Some(g), None) =>
        // 			(None, Some(b)) =>
        // 			(None, None) => panic!("Integrator builder for Newtonian mode requires at least one force.")
        // 		}
        // 	},
        // 	Some("Langevin") => {

        // 	},
        // 	None => panic!("Integrator builder requires mode. ex) Newton, Langevin")
        // }
    }
}

pub struct RKIntegrator<'a, V, S, T>(Box<dyn RKIntegratorTrait<'a, V, S, T>>);

impl<'a, V, S, T> RKIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    pub fn iterate(&mut self) {
        self.0.iterate()
    }

    fn from_builder<I>(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self
    where
        I: RKIntegratorTrait<'a, V, S, T>,
    {
        Self(Box::new(I::from_builder(builder)))
    }
}

pub struct NewtonGlobalIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: T,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_force: Vec<V>,
    temp_state: Vec<S>,
    global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> NewtonGlobalIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        match builder.mode {
            Some("Newton") => {}
            _ => unreachable!(),
        }
        let tableau = match builder.method {
            Some(m) => m.tableau().clone(),
            None => panic!("Integrator builder requires Runge-Kutta method information"),
        };
        let num_mem = tableau.len() + 1;

        let time_iter = match builder.time_iter {
            Some(t) => t.clone(),
            None => panic!("Integrator Builder requires time iterator."),
        };

        let states = match builder.states {
            Some(s) => s.clone(),
            None => panic!("Integrator Builder requires vector of states to iterate."),
        };
        let temp_state = states.clone();
        let num_states = states.len();
        if num_states == 0 {
            panic!("Integrator requires at least one state");
        }

        let zeros = states[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
        let memory = vec![vec![(zeros.clone(), zeros.clone()); num_states]; num_mem];
        let temp_force = vec![zeros.clone(); num_states];

        let global_force = match builder.global_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let boundary_cond = match builder.boundary_cond {
            Some(b) => b,
            None => panic!("Integrator builder requires boundary conditions."),
        };

        Self {
            tableau,
            time_iter,
            states,
            memory,
            temp_force,
            temp_state,
            global_force,
            boundary_cond,
        }
    }
}

impl<'a, V, S, T> RKIntegratorTrait for NewtonGlobalIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    fn force(&mut self) {
        unimplemented!();
    }
}

pub struct NewtonBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: T,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_force: Vec<V>,
    temp_state: Vec<S>,
    bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> NewtonBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        match builder.mode {
            Some("Newton") => {}
            _ => unreachable!(),
        }
        let tableau = match builder.method {
            Some(m) => m.tableau().clone(),
            None => panic!("Integrator builder requires Runge-Kutta method information"),
        };
        let num_mem = tableau.len() + 1;

        let time_iter = match builder.time_iter {
            Some(t) => t.clone(),
            None => panic!("Integrator Builder requires time iterator."),
        };

        let states = match builder.states {
            Some(s) => s.clone(),
            None => panic!("Integrator Builder requires vector of states to iterate."),
        };
        let temp_state = states.clone();
        let num_states = states.len();
        if num_states == 0 {
            panic!("Integrator requires at least one state");
        }

        let zeros = states[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
        let memory = vec![vec![(zeros.clone(), zeros.clone()); num_states]; num_mem];
        let temp_force = vec![zeros.clone(); num_states];

        let bimolecular_force = match builder.bimolecular_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let boundary_cond = match builder.boundary_cond {
            Some(b) => b,
            None => panic!("Integrator builder requires boundary conditions."),
        };

        Self {
            tableau,
            time_iter,
            states,
            memory,
            temp_force,
            temp_state,
            bimolecular_force,
            boundary_cond,
        }
    }
}

impl<'a, V, S, T> RKIntegratorTrait for NewtonBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn force(&mut self) {
        unimplemented!();
    }

    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        unimplemented!();
    }
}

pub struct NewtonBothIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: T,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_force: Vec<V>,
    temp_state: Vec<S>,
    global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> NewtonBothIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        match builder.mode {
            Some("Newton") => {}
            _ => unreachable!(),
        }
        let tableau = match builder.method {
            Some(m) => m.tableau().clone(),
            None => panic!("Integrator builder requires Runge-Kutta method information"),
        };
        let num_mem = tableau.len() + 1;

        let time_iter = match builder.time_iter {
            Some(t) => t.clone(),
            None => panic!("Integrator Builder requires time iterator."),
        };

        let states = match builder.states {
            Some(s) => s.clone(),
            None => panic!("Integrator Builder requires vector of states to iterate."),
        };
        let temp_state = states.clone();
        let num_states = states.len();
        if num_states == 0 {
            panic!("Integrator requires at least one state");
        }

        let zeros = states[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
        let memory = vec![vec![(zeros.clone(), zeros.clone()); num_states]; num_mem];
        let temp_force = vec![zeros.clone(); num_states];

        let global_force = match builder.global_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let bimolecular_force = match builder.bimolecular_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let boundary_cond = match builder.boundary_cond {
            Some(b) => b,
            None => panic!("Integrator builder requires boundary conditions."),
        };

        Self {
            tableau,
            time_iter,
            states,
            memory,
            temp_force,
            temp_state,
            global_force,
            bimolecular_force,
            boundary_cond,
        }
    }
}

impl<'a, V, S, T> RKIntegratorTrait for NewtonBothIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    fn force(&mut self) {
        unimplemented!();
    }
}

pub struct LangevinIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: T,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_force: Vec<V>,
    temp_state: Vec<S>,
    random_force: Box<dyn RandomForce<'a, S, Force = V>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> LangevinIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        match builder.mode {
            Some("Langevin") => {}
            _ => unreachable!(),
        }
        let tableau = match builder.method {
            Some(m) => m.tableau().clone(),
            None => panic!("Integrator builder requires Runge-Kutta method information"),
        };
        let num_mem = tableau.len() + 1;

        let time_iter = match builder.time_iter {
            Some(t) => t.clone(),
            None => panic!("Integrator Builder requires time iterator."),
        };

        let states = match builder.states {
            Some(s) => s.clone(),
            None => panic!("Integrator Builder requires vector of states to iterate."),
        };
        let temp_state = states.clone();
        let num_states = states.len();
        if num_states == 0 {
            panic!("Integrator requires at least one state");
        }

        let zeros = states[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
        let memory = vec![vec![(zeros.clone(), zeros.clone()); num_states]; num_mem];
        let temp_force = vec![zeros.clone(); num_states];

        let random_force = match builder.random_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let boundary_cond = match builder.boundary_cond {
            Some(b) => b,
            None => panic!("Integrator builder requires boundary conditions."),
        };

        Self {
            tableau,
            time_iter,
            states,
            memory,
            temp_force,
            temp_state,
            random_force,
            boundary_cond,
        }
    }
}

impl<'a, V, S, T> RKIntegratorTrait for LangevinIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    fn force(&mut self) {
        unimplemented!();
    }
}

pub struct LangevinGlobalIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: T,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_force: Vec<V>,
    temp_state: Vec<S>,
    global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    random_force: Box<dyn RandomForce<'a, S, Force = V>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> LangevinGlobalIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        match builder.mode {
            Some("Langevin") => {}
            _ => unreachable!(),
        }
        let tableau = match builder.method {
            Some(m) => m.tableau().clone(),
            None => panic!("Integrator builder requires Runge-Kutta method information"),
        };
        let num_mem = tableau.len() + 1;

        let time_iter = match builder.time_iter {
            Some(t) => t.clone(),
            None => panic!("Integrator Builder requires time iterator."),
        };

        let states = match builder.states {
            Some(s) => s.clone(),
            None => panic!("Integrator Builder requires vector of states to iterate."),
        };
        let temp_state = states.clone();
        let num_states = states.len();
        if num_states == 0 {
            panic!("Integrator requires at least one state");
        }

        let zeros = states[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
        let memory = vec![vec![(zeros.clone(), zeros.clone()); num_states]; num_mem];
        let temp_force = vec![zeros.clone(); num_states];

        let global_force = match builder.global_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let random_force = match builder.random_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let boundary_cond = match builder.boundary_cond {
            Some(b) => b,
            None => panic!("Integrator builder requires boundary conditions."),
        };

        Self {
            tableau,
            time_iter,
            states,
            memory,
            temp_force,
            temp_state,
            global_force,
            random_force,
            boundary_cond,
        }
    }
}

impl<'a, V, S, T> RKIntegratorTrait for LangevinGlobalIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    fn force(&mut self) {
        unimplemented!();
    }
}

pub struct LangevinBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: T,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_force: Vec<V>,
    temp_state: Vec<S>,
    bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    random_force: Box<dyn RandomForce<'a, S, Force = V>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> LangevinBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        match builder.mode {
            Some("Langevin") => {}
            _ => unreachable!(),
        }
        let tableau = match builder.method {
            Some(m) => m.tableau().clone(),
            None => panic!("Integrator builder requires Runge-Kutta method information"),
        };
        let num_mem = tableau.len() + 1;

        let time_iter = match builder.time_iter {
            Some(t) => t.clone(),
            None => panic!("Integrator Builder requires time iterator."),
        };

        let states = match builder.states {
            Some(s) => s.clone(),
            None => panic!("Integrator Builder requires vector of states to iterate."),
        };
        let temp_state = states.clone();
        let num_states = states.len();
        if num_states == 0 {
            panic!("Integrator requires at least one state");
        }

        let zeros = states[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
        let memory = vec![vec![(zeros.clone(), zeros.clone()); num_states]; num_mem];
        let temp_force = vec![zeros.clone(); num_states];

        let bimolecular_force = match builder.bimolecular_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let random_force = match builder.random_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let boundary_cond = match builder.boundary_cond {
            Some(b) => b,
            None => panic!("Integrator builder requires boundary conditions."),
        };

        Self {
            tableau,
            time_iter,
            states,
            memory,
            temp_force,
            temp_state,
            bimolecular_force,
            random_force,
            boundary_cond,
        }
    }
}

impl<'a, V, S, T> RKIntegratorTrait for LangevinBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    fn force(&mut self) {
        unimplemented!();
    }
}

pub struct LangevinBothIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: T,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_force: Vec<V>,
    temp_state: Vec<S>,
    global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    random_force: Box<dyn RandomForce<'a, S, Force = V>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> LangevinBothIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V>,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Clone,
    <V as Vector>::Item: Scalar,
{
    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
        match builder.mode {
            Some("Langevin") => {}
            _ => unreachable!(),
        }
        let tableau = match builder.method {
            Some(m) => m.tableau().clone(),
            None => panic!("Integrator builder requires Runge-Kutta method information"),
        };
        let num_mem = tableau.len() + 1;

        let time_iter = match builder.time_iter {
            Some(t) => t.clone(),
            None => panic!("Integrator Builder requires time iterator."),
        };

        let states = match builder.states {
            Some(s) => s.clone(),
            None => panic!("Integrator Builder requires vector of states to iterate."),
        };
        let temp_state = states.clone();
        let num_states = states.len();
        if num_states == 0 {
            panic!("Integrator requires at least one state");
        }

        let zeros = states[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
        let memory = vec![vec![(zeros.clone(), zeros.clone()); num_states]; num_mem];
        let temp_force = vec![zeros.clone(); num_states];

        let global_force = match builder.global_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let bimolecular_force = match builder.bimolecular_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let random_force = match builder.random_force {
            Some(f) => f,
            None => unreachable!(),
        };

        let boundary_cond = match builder.boundary_cond {
            Some(b) => b,
            None => panic!("Integrator builder requires boundary conditions."),
        };

        Self {
            tableau,
            time_iter,
            states,
            memory,
            temp_force,
            temp_state,
            global_force,
            bimolecular_force,
            random_force,
            boundary_cond,
        }
    }
}

impl<'a, V, S, T> RKIntegratorTrait for LangevinBothIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)>,
    T: TimeIterator<<V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    fn force(&mut self) {
        unimplemented!();
    }
}

trait RKIntegratorTrait<'a, V, S, T> {
    fn iterate(&mut self) {
        unimplemented!();
    }

    fn force(&mut self);

    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self
    where
        V: Vector,
        S: State<Position = V, Movement = (V, V)>,
        T: TimeIterator<<V as Vector>::Item>,
        <V as Vector>::Item: Scalar,
        Self: Sized;
}
