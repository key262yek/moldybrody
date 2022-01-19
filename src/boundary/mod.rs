


pub trait Boundary<V>{
    fn check_inclusion(&self, pos : &V) -> bool;
    fn normal_at(&self, pos : &V) -> Option<V>;
}


pub mod plane;
pub mod sphere;
