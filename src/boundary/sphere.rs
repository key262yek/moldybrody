

use crate::vector::Vector;

pub struct Sphere<V : Vector>{
    center : V,
    radius : <V as Vector>::Item,
}
