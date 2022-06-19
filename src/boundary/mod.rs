



use std::cell::RefCell;
use std::rc::Rc;
use crate::state::{State, HasVelocity};
use crate::vector::Vector;
use std::fmt::Debug;

use rand_pcg::Pcg64;


// #[derive(Debug)]
// pub struct System<V>{
//     #[allow(dead_code)]
//     bcs : Vec<(Box<dyn Boundary<V>>, BCspecies)>,
// }

// impl<V, T> System<V>
//     where V : Vector<Item = T>{
//     pub fn check_bc<'a, S, M>(&'a self, state : &'a S, movement : &'a mut M) -> AfterMove
//         where S : State<Movement = M, Position = V>,
//               M : Vector<Item = T>,
//               &'a V : Add<&'a V, Output = V>{

//         for (boundary, bc) in &self.bcs{
//             if let Some(intersect) = boundary.find_intersect(state.pos(), state.disp(movement)){
//                 unimplemented!();

//             }
//         }
//         return AfterMove::Survived;
//     }
// }

pub trait Boundary<V> : Debug{
    fn check_inclusion(&self, pos : &V) -> bool;

    // fn check_bc(&self, bc : BoundaryCondition, pos : &V, movement : &mut V) -> AfterMove{
    //     unimplemented!();
    // }

    // Normal vector to inside
    fn normal_at(&self, pos : &V) -> Option<V>;

    // function believe position is at boundary
    fn normal_at_unsafe(&self, pos : &V) -> V;

    fn find_intersect(&self, pos : &V, movement : &V) -> Option<V>;

    // function believe pos is in system,
    fn find_intersect_unsafe(&self, pos : &V, movement : &V) -> Option<V>;
}

pub trait FloatBoundary<V>
    where V : Vector{
    fn check_inclusion(&self, pos : &V) -> bool;

    // Normal vector to inside
    fn normal_at(&self, pos : &V) -> Option<V>;

    // function believe position is at boundary
    fn normal_at_unsafe(&self, pos : &V) -> V;

    fn ratio_to_intersect(&self, pos : &V, movement : &V) -> Option<<V as Vector>::Item>;

    fn ratio_to_intersect_unsafe(&self, pos : &V, movement : &V) -> Option<<V as Vector>::Item>;

    fn find_intersect(&self, pos : &V, movement : &V) -> Option<V>;

    // function believe pos is in system,
    fn find_intersect_unsafe(&self, pos : &V, movement : &V) -> Option<V>;
}

pub trait IntBoundary<V>
    where V : Vector{
    fn check_inclusion(&self, pos : &V) -> bool;

    // Normal vector to inside
    fn normal_at(&self, pos : &V) -> Option<V>;

    // function believe position is at boundary
    fn normal_at_unsafe(&self, pos : &V) -> V;

    fn find_intersect(&self, pos : &V, movement : &V) -> Option<V>;

    // function believe pos is in system,
    fn find_intersect_unsafe(&self, pos : &V, movement : &V) -> Option<V>;
}

pub trait NonPeriodic {}

pub trait Periodic<V : Vector>{
    fn find_pair(&self, pos : &V) -> V;

    fn find_pair_mut(&self, pos : &mut V);
}

pub trait NonPeriodicBoundaryCondition<V>
    where V : Vector{
    fn check_bc<'a, S>(&self, state : &'a S, movement : &'a mut (V, V)) -> AfterMove
        where S : State<Movement = (V, V), Position = V> + HasVelocity;

    fn check_bc_overdamped<'a, S>(&self, state : &'a S, movement : &'a mut V) -> AfterMove
        where S : State<Movement = V, Position = V>;
}

pub trait PeriodicBoundaryCondition<V>
    where V : Vector{
    fn check_bc<'a, S>(&self, state : &'a S, movement : &'a mut (V, V)) -> AfterMove
        where S : State<Movement = (V, V), Position = V> + HasVelocity;

    fn check_bc_overdamped<'a, S>(&self, state : &'a S, movement : &'a mut V) -> AfterMove
        where S : State<Movement = V, Position = V>;
}

#[derive(Debug, Clone)]
pub enum BCspecies<T>{
    Periodic,
    Absorbing,
    PartiallyAbsorbing(T, Rc<RefCell<Pcg64>>),
    Reflective,
    InelasticReflection(T),
    DiffusiveReflection,
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum AfterMove{
    Survived,
    Dead,
}


pub mod plane;
pub mod sphere;
pub mod boundarycondition;
