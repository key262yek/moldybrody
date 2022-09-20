use crate::vector::Vector;
use crate::vector::basic::Map;

impl<V: Vector, W: Vector> Vector for (V, W) {
    type Item = V::Item;

    fn dim(&self) -> usize {
        self.0.dim() + self.1.dim()
    }
}





#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::Cartessian;

    #[test]
    fn test_vector_trait_for_pair(){
        let v = Cartessian::new([2; 2]);
        let pair = (v.clone(), v.clone());

        assert_eq!(pair.dim(), 4);
    }
}