use crate::state::State;
use crate::vector::Vector;
use rand_pcg::Pcg64;
use std::fmt::Debug;
use std::sync::Arc;
use std::sync::RwLock;

pub trait FloatBoundary<V>: Debug
where
    V: Vector,
{
    fn check_inclusion(&self, pos: &V) -> bool;

    // Normal vector to inside
    fn normal_at(&self, pos: &V) -> Option<V>;

    // function believe position is at boundary
    fn normal_at_unsafe(&self, pos: &V) -> V;

    fn ratio_to_intersect(&self, pos: &V, movement: &V) -> Option<<V as Vector>::Item>;

    fn ratio_to_intersect_unsafe(&self, pos: &V, movement: &V) -> Option<<V as Vector>::Item>;

    fn find_intersect(&self, pos: &V, movement: &V) -> Option<V>;

    // function believe pos is in system,
    fn find_intersect_unsafe(&self, pos: &V, movement: &V) -> Option<V>;
}

pub trait IntBoundary<V>: Debug
where
    V: Vector,
{
    fn check_inclusion(&self, pos: &V) -> bool;

    // Normal vector to inside
    fn normal_at(&self, pos: &V) -> Option<V>;

    // function believe position is at boundary
    fn normal_at_unsafe(&self, pos: &V) -> V;

    fn find_intersect(&self, pos: &V, movement: &V) -> Option<V>;

    // function believe pos is in system,
    fn find_intersect_unsafe(&self, pos: &V, movement: &V) -> Option<V>;
}

pub trait NonPeriodic {}

pub trait Periodic<V: Vector> {
    fn find_pair(&self, pos: &V) -> V;

    fn find_pair_mut(&self, pos: &mut V);
}

#[derive(Debug, Clone)]
pub enum BCspecies<T> {
    Periodic,
    Absorbing,
    PartiallyAbsorbing(T, Arc<RwLock<Pcg64>>),
    Reflective,
    InelasticReflection(T),
    DiffusiveReflection,
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum AfterMove {
    Survived,
    Dead,
}

pub trait BoundaryCondition<S, V> : Debug
where
    S: State<Position = V, Movement = (V, V)>,
    V: Vector,
{
    fn check_bc<'a>(&self, state: &'a S, movement: &'a mut (V, V)) -> AfterMove;
}

pub trait OverdampedBoundaryCondition<S, V>
where
    S: State<Position = V, Movement = V>,
    V: Vector,
{
    fn check_bc_overdamped<'a>(&self, state: &'a S, movement: &'a mut V) -> AfterMove;
}

#[derive(Debug)]
pub struct System<S, V>
where
    S: State<Position = V, Movement = (V, V)>,
    V: Vector,
{
    pub bcs: Vec<Box<dyn BoundaryCondition<S, V>>>,
}

impl<S, V> System<S, V>
where
    S: State<Position = V, Movement = (V, V)>,
    V: Vector,
{
    pub fn check_bc<'a>(&self, state: &'a S, movement: &'a mut (V, V)) -> AfterMove {
        for bc in &self.bcs {
            match bc.check_bc(state, movement) {
                AfterMove::Survived => {
                    continue;
                }
                AfterMove::Dead => return AfterMove::Dead,
            }
        }

        return AfterMove::Survived;
    }
}

pub struct OverdampedSystem<S, V>
where
    S: State<Position = V, Movement = V>,
    V: Vector,
{
    pub bcs: Vec<Box<dyn OverdampedBoundaryCondition<S, V>>>,
}

impl<S, V> OverdampedSystem<S, V>
where
    S: State<Position = V, Movement = V>,
    V: Vector,
{
    pub fn check_bc_overdamped<'a>(&self, state: &'a S, movement: &'a mut V) -> AfterMove {
        for bc in &self.bcs {
            match bc.check_bc_overdamped(state, movement) {
                AfterMove::Survived => {
                    continue;
                }
                AfterMove::Dead => return AfterMove::Dead,
            }
        }

        return AfterMove::Survived;
    }
}

pub mod boundary_serde;
pub mod boundarycondition;
pub mod initialization;
pub mod plane;
pub mod sphere;
