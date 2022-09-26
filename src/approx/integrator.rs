use crate::approx::tableau::ApproxMethod;
use crate::approx::tableau::ButcherTableau;
use crate::approx::{TimeDiffIterator, TimeIterator};
use crate::boundary::AfterMove;
use crate::boundary::System;
use crate::force::{Bimolecular, Global, RandomForce};
use crate::state::State;
use crate::vector::basic::Map;
use crate::vector::Scalar;
use crate::vector::Vector;
use num_traits::Zero;
use std::fmt::Debug;
use std::ops::Mul;
use crate::state::{Mass, HasVelocity};
use std::sync::{Mutex, Arc};
use std::cell::{Cell, RefCell};
use std::rc::Rc;

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
        V: Clone
            + Mul<<V as Vector>::Item, Output = V>
            + Debug
            + Map<Item = <V as Vector>::Item>
            + 'static,
        S: Clone + Debug + 'static,
        Vec<S>: Clone,
        T: Clone + Debug + 'static,
    {
        // match self.mode {
        //     Some("Newton") => {
        //         match (
        //             self.global_force.is_some(),
        //             self.bimolecular_force.is_some(),
        //         ) {
        //             (true, true) => {
        //                 RKIntegrator::from_builder::<NewtonBothIntegrator<V, S, T>>(self)
        //             }
        //             (true, false) => {
        //                 RKIntegrator::from_builder::<NewtonGlobalIntegrator<V, S, T>>(self)
        //             }
        //             (false, true) => {
        //                 RKIntegrator::from_builder::<NewtonBimoleculeIntegrator<V, S, T>>(self)
        //             }
        //             (false, false) => {
        //                 panic!("Integrator builder for Newtonian mode requires at least one force.")
        //             }
        //         }
        //     }
        //     Some("Langevin") => {
        //         match (
        //             self.global_force.is_some(),
        //             self.bimolecular_force.is_some(),
        //             self.random_force.is_some(),
        //         ) {
        //             (true, true, true) => {
        //                 RKIntegrator::from_builder::<LangevinBothIntegrator<V, S, T>>(self)
        //             }
        //             (true, false, true) => {
        //                 RKIntegrator::from_builder::<LangevinGlobalIntegrator<V, S, T>>(self)
        //             }
        //             (false, true, true) => {
        //                 RKIntegrator::from_builder::<LangevinBimoleculeIntegrator<V, S, T>>(self)
        //             }
        //             (false, false, true) => {
        //                 RKIntegrator::from_builder::<LangevinIntegrator<V, S, T>>(self)
        //             }
        //             (_, _, false) => {
        //                 panic!("Integrator builder for Langevin mode requires random force")
        //             }
        //         }
        //     }
        //     Some(_) => panic!("Invalid mode input. only Newton, Langevin mode are provided"),
        //     None => panic!("Integrator builder requires mode. ex) Newton, Langevin"),
        // }
        unimplemented!();
    }
}

trait RKIntegratorTrait<'a, V, S, T>: Debug
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar{

    fn iterate(&mut self);

    fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self
        where Self : Sized;
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

macro_rules! impl_from_builder {
    ($mode : tt, $force : tt, "code") => {
        fn from_builder(builder: RKIntegratorBuilder<'a, V, S, T>) -> Self{
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

            let given = match builder.states {
                Some(s) => s.clone(),
                None => panic!("Integrator Builder requires vector of states to iterate."),
            };
            let states = Rc::new(given.iter().map(|state| RefCell::new(state.clone())).collect::<Vec<RefCell<S>>>());
            let temp_states = Rc::new(given.iter().map(|state| RefCell::new(state.clone())).collect::<Vec<RefCell<S>>>());
            let num_states = given.len();
            if num_states == 0 {
                panic!("Integrator requires at least one state");
            }

            let zeros = given[0].pos().clone() * <<V as Vector>::Item as Zero>::zero();
            let temp_forces = Rc::new((0..num_states).map(|_| RefCell::new(zeros.clone())).collect::<Vec<RefCell<V>>>());
            let memories = Rc::new((0..num_mem).map(|_| {
                            (0..num_states).map(|_| RefCell::new((zeros.clone(), zeros.clone()))).collect::<Vec<RefCell<(V,V)>>>()
                        }).collect::<Vec<Vec<RefCell<(V, V)>>>>());
            let aftermoves = Rc::new((0..num_states).map(|_| Cell::new(AfterMove::Survived)).collect::<Vec<Cell<AfterMove>>>());

            let boundary_cond = match builder.boundary_cond {
                Some(b) => b,
                None => panic!("Integrator builder requires boundary conditions."),
            };

            impl_from_builder!(
                $mode,
                $force,
                builder,
                tableau,
                time_iter,
                states,
                memories,
                temp_forces,
                temp_states,
                boundary_cond,
                aftermoves,
                Self
            );
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
    ("Newton", "Global", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memories : ident, $temp_forces : ident, $temp_states : ident, $boundary_cond : ident, $aftermoves : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memories,
            $temp_forces,
            $temp_states,
            $aftermoves,
            global_force: match $builder.global_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Newton", "Bimolecule", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memories : ident, $temp_forces : ident, $temp_states : ident, $boundary_cond : ident, $aftermoves : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memories,
            $temp_forces,
            $temp_states,
            $aftermoves,
            bimolecular_force: match $builder.bimolecular_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Newton", "Both", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memories : ident, $temp_forces : ident, $temp_states : ident, $boundary_cond : ident, $aftermoves : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memories,
            $temp_forces,
            $temp_states,
            $aftermoves,
            global_force: match $builder.global_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            bimolecular_force: match $builder.bimolecular_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memories : ident, $temp_forces : ident, $temp_states : ident, $boundary_cond : ident, $aftermoves : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memories,
            $temp_forces,
            $temp_states,
            $aftermoves,
            random_force: match $builder.random_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "Global", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memories : ident, $temp_forces : ident, $temp_states : ident, $boundary_cond : ident, $aftermoves : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memories,
            $temp_forces,
            $temp_states,
            $aftermoves,
            global_force: match $builder.global_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            random_force: match $builder.random_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "Bimolecule", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memories : ident, $temp_forces : ident, $temp_states : ident, $boundary_cond : ident, $aftermoves : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memories,
            $temp_forces,
            $temp_states,
            $aftermoves,
            bimolecular_force: match $builder.bimolecular_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            random_force: match $builder.random_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            $boundary_cond,
        }
    };
    ("Langevin", "Both", $builder : ident, $tableau : ident, $time_iter : ident, $states : ident, $memories : ident, $temp_forces : ident, $temp_states : ident, $boundary_cond : ident, $aftermoves : ident, $self : ident) => {
        return $self {
            $tableau,
            $time_iter,
            $states,
            $memories,
            $temp_forces,
            $temp_states,
            $aftermoves,
            global_force: match $builder.global_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            bimolecular_force: match $builder.bimolecular_force {
                Some(f) => Rc::new(f),
                None => unreachable!(),
            },
            random_force: match $builder.random_force {
                Some(f) => Rc::new(f),
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
    pub states: Rc<Vec<RefCell<S>>>,
    memories: Rc<Vec<Vec<RefCell<(V, V)>>>>,
    temp_forces: Rc<Vec<RefCell<V>>>,
    temp_states: Rc<Vec<RefCell<S>>>,
    global_force: Rc<Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    boundary_cond: System<S, V>,
    aftermoves: Rc<Vec<Cell<AfterMove>>>,
}

impl<'a, V, S, T> NewtonGlobalIntegrator<'a, V, S, T>
where
    V: Vector + Clone + Mul<<V as Vector>::Item, Output = V> + Debug + Map<Item = <V as Vector>::Item>,
    S: State<Position = V, Movement = (V, V)> + Mass<<V as Vector>::Item> + HasVelocity + Clone + Debug,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
    <V as Vector>::Item: Scalar,
{
    fn clone_states(&self){
        self.temp_states.iter().zip(self.states.iter()).for_each(|(t, s)| {
            let mut temp = t.borrow_mut();
            let state = s.borrow();
            temp.clone_from(&*state);
        });
    }

    fn compute_temp_states(&self, i : usize, h : <V as Vector>::Item, ci : <V as Vector>::Item, ai : &Vec<<V as Vector>::Item>){
        self.clone_states();
        self.memories.iter().take(i).zip(ai.iter().take(i)).for_each(|(memvec, aij)| {
            self.temp_states.iter().zip(memvec).for_each(|(t, m)| {
                let mut temp = t.borrow_mut();
                let memory = m.borrow();
                
                temp.renew_with_constant(&memory, *aij * h);
            })
        })
    }

    fn compute_force(&self){
        self.temp_states.iter().zip(self.temp_forces.iter()).for_each(|(s, f)| {
            let state = s.borrow();
            let mut force = f.borrow_mut();
            let forcefield = self.global_force.clone();

            let mass = state.mass();
            forcefield.force_to(&state, &mut force);
        });
    }
}


impl<'a, V, S, T> RKIntegratorTrait<'a, V, S, T> for NewtonGlobalIntegrator<'a, V, S, T>
where
    V: Vector + Clone + Mul<<V as Vector>::Item, Output = V> + Debug + Map<Item = <V as Vector>::Item>,
    S: State<Position = V, Movement = (V, V)> + Mass<<V as Vector>::Item> + HasVelocity + Clone + Debug,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item> + Clone + Debug,
    <V as Vector>::Item: Scalar,
{
    fn iterate(&mut self) {
        if let Some((t, h)) = self.time_iter.next() {
            let num_state = self.states.len();

            for (i, (ci, ai)) in self.tableau.into_iter().enumerate(){
                self.compute_temp_states(i, h, ci, ai);
                
                // for k in 0..num_state{
                //     let mass = self.temp_states[k].mass();
                //     self.global_force.force_to(&self.temp_states[k], &mut self.temp_forces[k]);

                //     self.memories[i][k].0.clone_from(self.temp_states[k].vel());
                //     self.memories[i][k].1.zip_mut_with(&self.temp_forces[k], |v, f| *v = f / mass);
                // }
            }
            
        }
    }

    impl_from_builder!("Newton", "Global", "code");
}

#[derive(Debug)]
pub struct NewtonBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: TimeDiffIterator<T>,
    pub states: Rc<Vec<RefCell<S>>>,
    memories: Rc<Vec<Vec<RefCell<(V, V)>>>>,
    temp_forces: Rc<Vec<RefCell<V>>>,
    temp_states: Rc<Vec<RefCell<S>>>,
    bimolecular_force: Rc<Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    boundary_cond: System<S, V>,
    aftermoves: Rc<Vec<Cell<AfterMove>>>,
}


#[derive(Debug)]
pub struct NewtonBothIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: TimeDiffIterator<T>,
    pub states: Rc<Vec<RefCell<S>>>,
    memories: Rc<Vec<Vec<RefCell<(V, V)>>>>,
    temp_forces: Rc<Vec<RefCell<V>>>,
    temp_states: Rc<Vec<RefCell<S>>>,
    global_force: Rc<Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    bimolecular_force: Rc<Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    boundary_cond: System<S, V>,
    aftermoves: Rc<Vec<Cell<AfterMove>>>,
}


#[derive(Debug)]
pub struct LangevinIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: TimeDiffIterator<T>,
    pub states: Rc<Vec<RefCell<S>>>,
    memories: Rc<Vec<Vec<RefCell<(V, V)>>>>,
    temp_forces: Rc<Vec<RefCell<V>>>,
    temp_states: Rc<Vec<RefCell<S>>>,
    random_force: Rc<Box<dyn RandomForce<'a, S, Force = V>>>,
    boundary_cond: System<S, V>,
    aftermoves: Rc<Vec<Cell<AfterMove>>>,
}


#[derive(Debug)]
pub struct LangevinGlobalIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: TimeDiffIterator<T>,
    pub states: Rc<Vec<RefCell<S>>>,
    memories: Rc<Vec<Vec<RefCell<(V, V)>>>>,
    temp_forces: Rc<Vec<RefCell<V>>>,
    temp_states: Rc<Vec<RefCell<S>>>,
    global_force: Rc<Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    random_force: Rc<Box<dyn RandomForce<'a, S, Force = V>>>,
    boundary_cond: System<S, V>,
    aftermoves: Rc<Vec<Cell<AfterMove>>>,
}


#[derive(Debug)]
pub struct LangevinBimoleculeIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: TimeDiffIterator<T>,
    pub states: Rc<Vec<RefCell<S>>>,
    memories: Rc<Vec<Vec<RefCell<(V, V)>>>>,
    temp_forces: Rc<Vec<RefCell<V>>>,
    temp_states: Rc<Vec<RefCell<S>>>,
    bimolecular_force: Rc<Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    random_force: Rc<Box<dyn RandomForce<'a, S, Force = V>>>,
    boundary_cond: System<S, V>,
    aftermoves: Rc<Vec<Cell<AfterMove>>>,
}


#[derive(Debug)]
pub struct LangevinBothIntegrator<'a, V, S, T>
where
    V: Vector,
    S: State<Position = V, Movement = (V, V)> + Clone,
    T: TimeIterator<<V as Vector>::Item> + Iterator<Item = <V as Vector>::Item>,
    <V as Vector>::Item: Scalar,
{
    tableau: ButcherTableau<<V as Vector>::Item>,
    time_iter: TimeDiffIterator<T>,
    pub states: Rc<Vec<RefCell<S>>>,
    memories: Rc<Vec<Vec<RefCell<(V, V)>>>>,
    temp_forces: Rc<Vec<RefCell<V>>>,
    temp_states: Rc<Vec<RefCell<S>>>,
    global_force: Rc<Box<dyn Global<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    bimolecular_force: Rc<Box<dyn Bimolecular<'a, S, Force = V, Potential = <V as Vector>::Item>>>,
    random_force: Rc<Box<dyn RandomForce<'a, S, Force = V>>>,
    boundary_cond: System<S, V>,
    aftermoves: Rc<Vec<Cell<AfterMove>>>,
}


#[cfg(test)]
mod test {
    use super::*;
    use crate::approx::tableau::NewtonEulerMethod;
    use crate::approx::time::ConstStep;
    use crate::boundary::System;
    use crate::force::global::ConstGravity;
    use crate::prelude::State;
    use crate::state::{HasVelocity, Mass};
    use crate::vector::Cartessian2D;
    use approx::assert_abs_diff_eq;

    #[derive(State, Clone, Debug)]
    struct TestState {
        mass: f64,
        pos: Cartessian2D<f64>,
        vel: Cartessian2D<f64>,
    }

    macro_rules! declare_builder{
        () => {
            NewtonGlobalIntegrator::from_builder(RKIntegratorBuilder::new()
                .mode("Newton")
                .method(Box::new(NewtonEulerMethod::<f64>::new()))
                .time_iter(ConstStep::<f64>::new(1e-3).unwrap())
                .states(vec![TestState {
                    mass: 1f64,
                    pos: Cartessian2D::<f64>::zeros(),
                    vel: Cartessian2D::<f64>::zeros(),
                }])
                .global_force(Box::new(ConstGravity::new(Cartessian2D::new([0f64, -10f64]))))
                .boundary_cond(System { bcs: vec![] }))
        }
    }

    #[test]
    fn test_clone_states(){
        let integrator = declare_builder!();

        {
            let mut state = integrator.states[0].lock().unwrap();
            let temp_states = integrator.temp_states[0].lock().unwrap();
            assert_abs_diff_eq!(state.pos(), temp_states.pos(), epsilon = 1e-3);

            state.pos[1] = -2.0;
        }
        integrator.clone_states();
        {
            let mut state = integrator.states[0].lock().unwrap();
            let temp_states = integrator.temp_states[0].lock().unwrap();
            
            assert_abs_diff_eq!(state.pos(), temp_states.pos(), epsilon = 1e-3);
        }
    }

    #[test]
    fn test_compute_temp_state(){
        let integrator = declare_builder!();
        let ci = 1.0;
        let ai = vec![1.0, 2.0];

        {
            for i in 0..2{
                let mut mem = integrator.memories[i][0].lock().unwrap();
                mem.0[0] = (i + 1) as f64;
            }
        }

        integrator.compute_temp_states(0, ci, &ai);
        {
            let temp = integrator.temp_states[0].lock().unwrap();
            assert_abs_diff_eq!(temp.pos[0], 0.0, epsilon = 1e-3);
        }

        integrator.compute_temp_states(1, ci, &ai);
        {
            let temp = integrator.temp_states[0].lock().unwrap();
            assert_abs_diff_eq!(temp.pos[0], 1.0, epsilon = 1e-3);
        }

        integrator.compute_temp_states(2, ci, &ai);
        {
            let temp = integrator.temp_states[0].lock().unwrap();
            assert_abs_diff_eq!(temp.pos[0], 5.0, epsilon = 1e-3);
        }

        
    }


    // #[test]
    // fn test_from_builder() {
    //     let mode = "Newton";
    //     let method = Box::new(NewtonEulerMethod::<f64>::new());
    //     let time_iter = ConstStep::<f64>::new(1e-3).unwrap();

    //     let states = vec![TestState {
    //         mass: 1f64,
    //         pos: Cartessian2D::<f64>::zeros(),
    //         vel: Cartessian2D::<f64>::zeros(),
    //     }];
    //     let force = Gravity::<f64>::new(1f64);

    //     let system: System<TestState, Cartessian2D<f64>> = System { bcs: vec![] };

    //     let builder = RKIntegratorBuilder::new()
    //         .mode(mode)
    //         .method(method)
    //         .time_iter(time_iter)
    //         .states(states)
    //         .bimolecular_force(Box::new(force))
    //         .boundary_cond(system);

    //     let integrator = NewtonBimoleculeIntegrator::from_builder(builder);

    //     assert_abs_diff_eq!(
    //         integrator.tableau,
    //         NewtonEulerMethod::<f64>::new().tableau()
    //     );
    //     assert_abs_diff_eq!(integrator.memories.len(), 2);
    //     assert_abs_diff_eq!(integrator.temp_forces.len(), 1);
    // }

    // #[test]
    // fn test_build() {
    //     let mode = "Newton";
    //     let method = Box::new(NewtonEulerMethod::<f64>::new());
    //     let time_iter = ConstStep::<f64>::new(1e-3).unwrap();

    //     let states = vec![TestState {
    //         mass: 1f64,
    //         pos: Cartessian2D::<f64>::zeros(),
    //         vel: Cartessian2D::<f64>::zeros(),
    //     }];
    //     let force = Gravity::<f64>::new(1f64);

    //     let system: System<TestState, Cartessian2D<f64>> = System { bcs: vec![] };

    //     let integrator = RKIntegratorBuilder::new()
    //         .mode(mode)
    //         .method(method)
    //         .time_iter(time_iter)
    //         .states(states)
    //         .bimolecular_force(Box::new(force))
    //         .boundary_cond(system)
    //         .build();
    // }
}
