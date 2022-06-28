use crate::vector::Dim;
use crate::vector::Scalar;
use crate::vector::Vector;

#[derive(Clone, Debug)]
pub enum InitShape<V, const N: usize>
where
    V: Vector + Dim<N> + Clone,
    <V as Vector>::Item: Scalar + Copy,
{
    UniformCube(V, <V as Vector>::Item),           // Radius
    UniformBox(V, [<V as Vector>::Item; N]),       // Center, Width's
    UniformSphere(V, <V as Vector>::Item),         // Center, Radius
    UniformSemiSphere(V, V, <V as Vector>::Item),  // Center, Norm, Radius
    GaussianSphere(V, <V as Vector>::Item),        // Center, Deviation
    GaussianSemiSphere(V, V, <V as Vector>::Item), // Center, Norm, Deviation
}
