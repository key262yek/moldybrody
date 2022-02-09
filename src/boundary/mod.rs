


use std::fmt::Debug;

#[derive(Debug)]
pub struct System<V>{
    bcs : Vec<(Box<dyn Boundary<V>>, BoundaryCondition)>,
}

// impl<V, T> System<V>
//     where V : Vector<Item = T>{
//     fn check_bc<'a, S, M>(&'a self, state : &'a S, movement : &'a mut M) -> AfterMove
//         where S : State<Movement = M, Position = V>,
//               M : Vector<Item = T>,
//               &'a V : Add<&'a V, Output = V>{

//         let destination = state.pos() + state.disp(movement);
//         for (boundary, bc) in self.bcs{
//             if !boundary.check_inclusion(&destination){

//             }
//         }
//     }
// }

pub trait Boundary<V> : Debug{
    fn check_inclusion(&self, pos : &V) -> bool;

    // Normal vector to inside
    fn normal_at(&self, pos : &V) -> Option<V>;

    // function believe position is at boundary
    fn normal_at_unsafe(&self, pos : &V) -> V;

    fn find_intersect(&self, pos : &V, movement : &V) -> Option<V>;

    // function believe pos is in system,
    fn find_intersect_unsafe(&self, pos : &V, movement : &V) -> Option<V>;
}

#[derive(Debug, Clone)]
pub enum BoundaryCondition{
    Periodic,
    Absorbing,
    PartiallyAbsorbing(f64),
    Reflective,
    InelasticReflection,
    DiffusiveReflection,
}

#[derive(Debug, Clone)]
pub enum AfterMove{
    Survived,
    Dead,
}


pub mod plane;
pub mod sphere;
pub mod boundarycondition;
