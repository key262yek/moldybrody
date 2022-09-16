use crate::approx::tableau::ApproxMethod;
use crate::approx::tableau::ButcherTableau;
use crate::approx::{TimeIterator, TimeDiffIterator};
use crate::boundary::System;
use crate::force::{Bimolecular, Global, RandomForce};
use crate::state::State;
use crate::vector::Scalar;
use crate::vector::Vector;
use crate::vector::basic::Map;
use num_traits::Zero;
use std::ops::Mul;
use std::fmt::Debug;

#[derive(Debug)]
pub struct RKIntegratorBuilder<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
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
impl<V, S, T> RKIntegratorBuilder<'static, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
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

    pub fn mode(mut self, mode: &'static str) -> Self {
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
        force: Box<dyn Global<'static, S, Force = V, Potential = <V as Vector>::Item>>,
    ) -> Self {
        self.global_force = Some(force);
        self
    }

    pub fn bimolecular_force(
        mut self,
        force: Box<dyn Bimolecular<'static, S, Force = V, Potential = <V as Vector>::Item>>,
    ) -> Self {
        self.bimolecular_force = Some(force);
        self
    }

    pub fn random_force(mut self, force: Box<dyn RandomForce<'static, S, Force = V>>) -> Self {
        self.random_force = Some(force);
        self
    }

    pub fn boundary_cond(mut self, sys: System<S, V>) -> Self {
        self.boundary_cond = Some(sys);
        self
    }

    pub fn build(self) -> RKIntegrator<'static, V, S, T> 
    where 
        V: std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug + Map<Item = <V as Vector>::Item> + 'static,
        S: std::clone::Clone + Debug + 'static,
        Vec<S>: Clone,
        T: Clone + Debug + 'static,
    {
        // match self.mode {
        // 	Some("Newton") => {
        // 		match (self.global_force.is_some(), self.bimolecular_force.is_some()) {
        // 			(true, true) => RKIntegrator::from_builder::<NewtonBothIntegrator<V,S,T>>(self),
        // 			(true, false) => RKIntegrator::from_builder::<NewtonGlobalIntegrator<V,S,T>>(self),
        // 			(false, true) => RKIntegrator::from_builder::<NewtonBimoleculeIntegrator<V,S,T>>(self),
        // 			(false, false) => panic!("Integrator builder for Newtonian mode requires at least one force.")
        // 		}
        // 	},
        // 	Some("Langevin") => {
        //         match (self.global_force.is_some(), self.bimolecular_force.is_some(), self.random_force.is_some()) {
        //             (true, true, true) => RKIntegrator::from_builder::<LangevinBothIntegrator<V,S,T>>(self),
        //             (true, false, true) => RKIntegrator::from_builder::<LangevinGlobalIntegrator<V,S,T>>(self),
        // 			(false, true, true) => RKIntegrator::from_builder::<LangevinBimoleculeIntegrator<V,S,T>>(self),
        //             (false, false, true) => RKIntegrator::from_builder::<LangevinIntegrator<V,S,T>>(self),
        //             (_, _, false) => panic!("Integrator builder for Langevin mode requires random force")
        // 		}
        // 	},
        //     Some(_) => panic!("Invalid mode input. only Newton, Langevin mode are provided"),
        // 	None => panic!("Integrator builder requires mode. ex) Newton, Langevin")
        // }
        unimplemented!();
    }
}

trait RKIntegratorTrait<'a, V, S, T> : Debug{
    fn iterate(&mut self) 
    where
        V: Vector,
        S: State<Position = V, Movement = (V, V)> + Clone,
        T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
        <V as Vector>::Item: Scalar,
    {
        if let Some((t, h)) = self.next_time(){
            self.reset_temp();
            // let length = self.tableau().len();

            // for i in 0..length{
            //     self.compute_force(i);
            //     self.compute_temp_state(i);
            // }
        }
    }

    fn compute_force(&'a mut self, idx : usize);
    
    fn compute_temp_state(&'a mut self, idx : usize);

    fn compute_next_state(&'a mut self);

    fn states_mut(&mut self) -> &mut Vec<S>;

    fn states(&self) -> &Vec<S>;

    fn reset_temp(&mut self);

    fn next_time(&mut self) -> Option<(<V as Vector>::Item, <V as Vector>::Item)>
        where V : Vector;

    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self
    where
        V: Vector,
        S: State<Position = V, Movement = (V, V)> + Clone,
        T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
        <V as Vector>::Item: Scalar,
        Self: Sized;
}

#[derive(Debug)]
pub struct RKIntegrator<'a, V, S, T>(Box<dyn RKIntegratorTrait<'a, V, S, T>>);

impl<'a, V, S, T> RKIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    pub fn iterate(&'a mut self) {
        self.0.iterate()
    }

    fn from_builder<I>(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self
    where
        I: RKIntegratorTrait<'a, V, S, T> + 'static,
    {
        Self(Box::new(I::from_builder(builder)))
    }
}

macro_rules! impl_from_builder{
    ($mode : tt, $force : tt, "code") => {
        fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self {
            impl_from_builder!($mode, builder);

            let tableau = match builder.method {
                Some(m) => m.tableau().clone(),
                None => panic!("Integrator builder requires Runge-Kutta method information"),
            };
            let num_mem = tableau.len() + 1;

            let time_iter = match builder.time_iter {
                Some(t) => t.into_diff(),
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

            let boundary_cond = match builder.boundary_cond {
                Some(b) => b,
                None => panic!("Integrator builder requires boundary conditions."),
            };

            impl_from_builder!($mode, $force, builder, tableau, time_iter, states, memory, temp_state, boundary_cond, Self);
        }
    };
    ("Newton", $builder : ident) => {
        match $builder.mode {
            Some("Newton") => {}
            _ => unreachable!(),
        }   
    };
    ("Langevin", $builder : ident) => {
        match $builder.mode {
            Some("Langevin") => {}
            _ => unreachable!(),
        }   
    };
    ("Newton", "Global", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memory : ident, $temp_state : ident, $boundary_cond : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memory,
            $temp_state,
            global_force : match $builder.global_force {
                Some(f) => f,
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Newton", "Bimolecule", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memory : ident, $temp_state : ident, $boundary_cond : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memory,
            $temp_state,
            bimolecular_force : match $builder.bimolecular_force {
                Some(f) => f,
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Newton", "Both", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memory : ident, $temp_state : ident, $boundary_cond : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memory,
            $temp_state,
            global_force : match $builder.global_force {
                Some(f) => f,
                None => unreachable!(),
            },
            bimolecular_force : match $builder.bimolecular_force {
                Some(f) => f,
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memory : ident, $temp_state : ident, $boundary_cond : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memory,
            $temp_state,
            random_force : match $builder.random_force {
                Some(f) => f,
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "Global", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memory : ident, $temp_state : ident, $boundary_cond : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memory,
            $temp_state,
            global_force : match $builder.global_force {
                Some(f) => f,
                None => unreachable!(),
            },
            random_force : match $builder.random_force {
                Some(f) => f,
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "Bimolecule", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memory : ident, $temp_state : ident, $boundary_cond : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memory,
            $temp_state,
            bimolecular_force : match $builder.bimolecular_force {
                Some(f) => f,
                None => unreachable!(),
            },
            random_force : match $builder.random_force {
                Some(f) => f,
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "Both", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memory : ident, $temp_state : ident, $boundary_cond : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memory,
            $temp_state,
            global_force : match $builder.global_force {
                Some(f) => f,
                None => unreachable!(),
            },
            bimolecular_force : match $builder.bimolecular_force {
                Some(f) => f,
                None => unreachable!(),
            },
            random_force : match $builder.random_force {
                Some(f) => f,
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
}


#[derive(Debug)]
pub struct NewtonGlobalIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: TimeDiffIterator<T>,
    pub states: Vec<S>,
    memory: Vec<Vec<(V, V)>>,
    temp_state: Vec<S>,
    global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
    boundary_cond: System<S, V>,
}

impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for NewtonGlobalIntegrator<'a, V, S, T>
where
    V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug,
    S: State<Position = V, Movement = (V, V)> + std::clone::Clone + Debug,
    Vec<S>: Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
    <V as Vector>::Item: Scalar,
{
    fn compute_force(&'a mut self, idx : usize) {
        self.temp_state.iter().zip(self.memory[idx].iter_mut()).for_each(|(state, f)| {
            self.global_force.force_to(state, &mut f.1);
        });
    }
    
    fn compute_temp_state(&'a mut self, idx : usize){
        unimplemented!();
    }

    fn compute_next_state(&'a mut self){
        unimplemented!();
    }

    fn states_mut(&mut self) -> &mut Vec<S>{
        &mut self.states
    }

    fn states(&self) -> &Vec<S>{
        &self.states
    }
    
    fn reset_temp(&mut self){
        self.temp_state.clone_from(&self.states);
    }

    fn next_time(&mut self) -> Option<(<V as Vector>::Item, <V as Vector>::Item)>{
        self.time_iter.next()
    }

    impl_from_builder!("Newton", "Global", "code");
}

// #[derive(Debug)]
// pub struct NewtonBimoleculeIntegrator<'a, V, S, T>
// where
//     V: Vector,
//     S: State<Position = V, Movement = (V, V)> + Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
//     <V as Vector>::Item: Scalar,
// {
//     tableau: ButcherTableau<<V as Vector>::Item>,
//     time_iter: TimeDiffIterator<T>,
//     pub states: Vec<S>,
//     memory: Vec<Vec<(V, V)>>,
//     temp_state: Vec<S>,
//     bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
//     boundary_cond: System<S, V>,
// }

// impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for NewtonBimoleculeIntegrator<'a, V, S, T>
// where
//     V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug + Map<Item = <V as Vector>::Item>,
//     S: State<Position = V, Movement = (V, V)> + std::clone::Clone + Debug,
//     Vec<S>: Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
//     <V as Vector>::Item: Scalar,
// {
//     fn compute_force(&'a mut self, idx : usize) {
//         unimplemented!();
//     }
    
//     fn compute_temp_state(&'a mut self, idx : usize){
//         unimplemented!();
//     }

//     fn compute_next_state(&'a mut self){
//         unimplemented!();
//     }

//     fn states_mut(&mut self) -> &mut Vec<S>{
//         &mut self.states
//     }

//     fn states(&self) -> &Vec<S>{
//         &self.states
//     }

//     fn time_iter(&mut self) -> &mut TimeDiffIterator<T>{
//         &mut self.time_iter
//     }

//     fn tableau(&self) -> &ButcherTableau<<V as Vector>::Item>{
//         &self.tableau
//     }

//     fn temp_state_mut(&mut self) -> &mut Vec<S>{
//         &mut self.temp_state
//     }

//     fn temp_state(&self) -> &Vec<S>{
//         &self.temp_state
//     }

//     impl_from_builder!("Newton", "Bimolecule", "code");
// }

// #[derive(Debug)]
// pub struct NewtonBothIntegrator<'a, V, S, T>
// where
//     V: Vector,
//     S: State<Position = V, Movement = (V, V)> + Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
//     <V as Vector>::Item: Scalar,
// {
//     tableau: ButcherTableau<<V as Vector>::Item>,
//     time_iter: TimeDiffIterator<T>,
//     pub states: Vec<S>,
//     memory: Vec<Vec<(V, V)>>,
//     temp_state: Vec<S>,
//     global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
//     bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
//     boundary_cond: System<S, V>,
// }

// impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for NewtonBothIntegrator<'a, V, S, T>
// where
//     V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug,
//     S: State<Position = V, Movement = (V, V)> + std::clone::Clone + Debug,
//     Vec<S>: Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
//     <V as Vector>::Item: Scalar,
// {
//     fn compute_force(&'a mut self, idx : usize) {
//         unimplemented!();
//     }
    
//     fn compute_temp_state(&'a mut self, idx : usize){
//         unimplemented!();
//     }

//     fn compute_next_state(&'a mut self){
//         unimplemented!();
//     }

//     fn states_mut(&mut self) -> &mut Vec<S>{
//         &mut self.states
//     }

//     fn states(&self) -> &Vec<S>{
//         &self.states
//     }

//     fn time_iter(&mut self) -> &mut TimeDiffIterator<T>{
//         &mut self.time_iter
//     }

//     fn tableau(&self) -> &ButcherTableau<<V as Vector>::Item>{
//         &self.tableau
//     }

//     fn temp_state_mut(&mut self) -> &mut Vec<S>{
//         &mut self.temp_state
//     }

//     fn temp_state(&self) -> &Vec<S>{
//         &self.temp_state
//     }

//     impl_from_builder!("Newton", "Both", "code");
// }

// #[derive(Debug)]
// pub struct LangevinIntegrator<'a, V, S, T>
// where
//     V: Vector,
//     S: State<Position = V, Movement = (V, V)> + Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
//     <V as Vector>::Item: Scalar,
// {
//     tableau: ButcherTableau<<V as Vector>::Item>,
//     time_iter: TimeDiffIterator<T>,
//     pub states: Vec<S>,
//     memory: Vec<Vec<(V, V)>>,
//     temp_state: Vec<S>,
//     random_force: Box<dyn RandomForce<'a, S, Force = V>>,
//     boundary_cond: System<S, V>,
// }

// impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for LangevinIntegrator<'a, V, S, T>
// where
//     V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug,
//     S: State<Position = V, Movement = (V, V)> + std::clone::Clone + Debug,
//     Vec<S>: Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
//     <V as Vector>::Item: Scalar,
// {
//     fn compute_force(&'a mut self, idx : usize) {
//         unimplemented!();
//     }
    
//     fn compute_temp_state(&'a mut self, idx : usize){
//         unimplemented!();
//     }

//     fn compute_next_state(&'a mut self){
//         unimplemented!();
//     }

//     fn states_mut(&mut self) -> &mut Vec<S>{
//         &mut self.states
//     }

//     fn states(&self) -> &Vec<S>{
//         &self.states
//     }

//     fn time_iter(&mut self) -> &mut TimeDiffIterator<T>{
//         &mut self.time_iter
//     }

//     fn tableau(&self) -> &ButcherTableau<<V as Vector>::Item>{
//         &self.tableau
//     }

//     fn temp_state_mut(&mut self) -> &mut Vec<S>{
//         &mut self.temp_state
//     }

//     fn temp_state(&self) -> &Vec<S>{
//         &self.temp_state
//     }

//     impl_from_builder!("Langevin", "", "code");
// }

// #[derive(Debug)]
// pub struct LangevinGlobalIntegrator<'a, V, S, T>
// where
//     V: Vector,
//     S: State<Position = V, Movement = (V, V)> + Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
//     <V as Vector>::Item: Scalar,
// {
//     tableau: ButcherTableau<<V as Vector>::Item>,
//     time_iter: TimeDiffIterator<T>,
//     pub states: Vec<S>,
//     memory: Vec<Vec<(V, V)>>,
//     temp_state: Vec<S>,
//     global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
//     random_force: Box<dyn RandomForce<'a, S, Force = V>>,
//     boundary_cond: System<S, V>,
// }

// impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for LangevinGlobalIntegrator<'a, V, S, T>
// where
//     V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug,
//     S: State<Position = V, Movement = (V, V)> + std::clone::Clone + Debug,
//     Vec<S>: Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
//     <V as Vector>::Item: Scalar,
// {
//     fn compute_force(&'a mut self, idx : usize) {
//         unimplemented!();
//     }
    
//     fn compute_temp_state(&'a mut self, idx : usize){
//         unimplemented!();
//     }

//     fn compute_next_state(&'a mut self){
//         unimplemented!();
//     }

//     fn states_mut(&mut self) -> &mut Vec<S>{
//         &mut self.states
//     }

//     fn states(&self) -> &Vec<S>{
//         &self.states
//     }

//     fn time_iter(&mut self) -> &mut TimeDiffIterator<T>{
//         &mut self.time_iter
//     }

//     fn tableau(&self) -> &ButcherTableau<<V as Vector>::Item>{
//         &self.tableau
//     }

//     fn temp_state_mut(&mut self) -> &mut Vec<S>{
//         &mut self.temp_state
//     }

//     fn temp_state(&self) -> &Vec<S>{
//         &self.temp_state
//     }

//     impl_from_builder!("Langevin", "Global", "code");
// }

// #[derive(Debug)]
// pub struct LangevinBimoleculeIntegrator<'a, V, S, T>
// where
//     V: Vector,
//     S: State<Position = V, Movement = (V, V)> + Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
//     <V as Vector>::Item: Scalar,
// {
//     tableau: ButcherTableau<<V as Vector>::Item>,
//     time_iter: TimeDiffIterator<T>,
//     pub states: Vec<S>,
//     memory: Vec<Vec<(V, V)>>,
//     temp_state: Vec<S>,
//     bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
//     random_force: Box<dyn RandomForce<'a, S, Force = V>>,
//     boundary_cond: System<S, V>,
// }

// impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for LangevinBimoleculeIntegrator<'a, V, S, T>
// where
//     V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug,
//     S: State<Position = V, Movement = (V, V)> + std::clone::Clone + Debug,
//     Vec<S>: Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
//     <V as Vector>::Item: Scalar,
// {
//     fn compute_force(&'a mut self, idx : usize) {
//         unimplemented!();
//     }
    
//     fn compute_temp_state(&'a mut self, idx : usize){
//         unimplemented!();
//     }

//     fn compute_next_state(&'a mut self){
//         unimplemented!();
//     }

//     fn states_mut(&mut self) -> &mut Vec<S>{
//         &mut self.states
//     }

//     fn states(&self) -> &Vec<S>{
//         &self.states
//     }

//     fn time_iter(&mut self) -> &mut TimeDiffIterator<T>{
//         &mut self.time_iter
//     }

//     fn tableau(&self) -> &ButcherTableau<<V as Vector>::Item>{
//         &self.tableau
//     }

//     fn temp_state_mut(&mut self) -> &mut Vec<S>{
//         &mut self.temp_state
//     }

//     fn temp_state(&self) -> &Vec<S>{
//         &self.temp_state
//     }

//     impl_from_builder!("Langevin", "Bimolecule", "code");
// }

// #[derive(Debug)]
// pub struct LangevinBothIntegrator<'a, V, S, T>
// where
//     V: Vector,
//     S: State<Position = V, Movement = (V, V)> + Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
//     <V as Vector>::Item: Scalar,
// {
//     tableau: ButcherTableau<<V as Vector>::Item>,
//     time_iter: TimeDiffIterator<T>,
//     pub states: Vec<S>,
//     memory: Vec<Vec<(V, V)>>,
//     temp_state: Vec<S>,
//     global_force: Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>,
//     bimolecular_force: Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>,
//     random_force: Box<dyn RandomForce<'a, S, Force = V>>,
//     boundary_cond: System<S, V>,
// }

// impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for LangevinBothIntegrator<'a, V, S, T>
// where
//     V: Vector + std::clone::Clone + Mul<<V as Vector>::Item, Output = V> + Debug,
//     S: State<Position = V, Movement = (V, V)> + std::clone::Clone + Debug,
//     Vec<S>: Clone,
//     T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
//     <V as Vector>::Item: Scalar,
// {
//     fn compute_force(&'a mut self, idx : usize) {
//         unimplemented!();
//     }
    
//     fn compute_temp_state(&'a mut self, idx : usize){
//         unimplemented!();
//     }

//     fn compute_next_state(&'a mut self){
//         unimplemented!();
//     }

//     fn states_mut(&mut self) -> &mut Vec<S>{
//         &mut self.states
//     }

//     fn states(&self) -> &Vec<S>{
//         &self.states
//     }

//     fn time_iter(&mut self) -> &mut TimeDiffIterator<T>{
//         &mut self.time_iter
//     }

//     fn tableau(&self) -> &ButcherTableau<<V as Vector>::Item>{
//         &self.tableau
//     }

//     fn temp_state_mut(&mut self) -> &mut Vec<S>{
//         &mut self.temp_state
//     }

//     fn temp_state(&self) -> &Vec<S>{
//         &self.temp_state
//     }

//     impl_from_builder!("Langevin", "Both", "code");
// }

// #[cfg(test)]
// mod test{
//     use super::*;
//     use crate::approx::tableau::NewtonEulerMethod;
//     use crate::approx::time::ConstStep;
//     use crate::prelude::State;
//     use crate::state::{Mass, HasVelocity};
//     use crate::vector::Cartessian2D;
//     use crate::force::bimolecular::Gravity;
//     use crate::boundary::System;
//     use approx::assert_abs_diff_eq;

//     #[test] 
//     fn test_from_builder(){
//         let mode = "Newton";
//         let method = Box::new(NewtonEulerMethod::<f64>::new());
//         let time_iter = ConstStep::<f64>::new(1e-3).unwrap();

//         #[derive(State, Clone, Debug)]
//         struct TestState{
//             mass : f64,
//             pos: Cartessian2D<f64>,
//             vel: Cartessian2D<f64>,
//         }

//         let states = vec![TestState{mass :1f64, pos : Cartessian2D::<f64>::zeros(), vel : Cartessian2D::<f64>::zeros()}];
//         let force = Gravity::<f64>::new(1f64);

//         let system : System<TestState, Cartessian2D<f64>> = System{bcs : vec![]};        

//         let builder = RKIntegratorBuilder::new()
//             .mode(mode)
//             .method(method)
//             .time_iter(time_iter)
//             .states(states)
//             .bimolecular_force(Box::new(force))
//             .boundary_cond(system);

//         let integrator = NewtonBimoleculeIntegrator::from_builder(builder);

//         assert_abs_diff_eq!(integrator.tableau, NewtonEulerMethod::<f64>::new().tableau());
//         assert_abs_diff_eq!(integrator.memory.len(), 2);
//         assert_abs_diff_eq!(integrator.temp_force.len(), 1);
//     }

//     #[test] 
//     fn test_build(){
//         let mode = "Newton";
//         let method = Box::new(NewtonEulerMethod::<f64>::new());
//         let time_iter = ConstStep::<f64>::new(1e-3).unwrap();

//         #[derive(State, Clone, Debug)]
//         struct TestState{
//             mass : f64,
//             pos: Cartessian2D<f64>,
//             vel: Cartessian2D<f64>,
//         }

//         let states = vec![TestState{mass :1f64, pos : Cartessian2D::<f64>::zeros(), vel : Cartessian2D::<f64>::zeros()}];
//         let force = Gravity::<f64>::new(1f64);

//         let system : System<TestState, Cartessian2D<f64>> = System{bcs : vec![]};        

//         let integrator = RKIntegratorBuilder::new()
//             .mode(mode)
//             .method(method)
//             .time_iter(time_iter)
//             .states(states)
//             .bimolecular_force(Box::new(force))
//             .boundary_cond(system)
//             .build();

        
//     }
// }