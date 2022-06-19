
use serde::{ser::{SerializeStruct, SerializeSeq}, Serialize, Serializer};
use super::{Cartessian, CartessianND, arithmetic::Scalar};
use std::slice::Iter;

struct Sequence<'a, T>(Iter<'a, T>);

impl<'a, T> Serialize for Sequence<'a, T>
    where T : Serialize{

    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: Serializer {
        let iter = &self.0;
        let mut seq = serializer.serialize_seq(Some(iter.len()))?;
        for elt in iter.clone(){
            seq.serialize_element(elt)?;
        }
        seq.end()
    }
}

// ========================================================

impl<T, const N : usize> Serialize for Cartessian<T, N>
    where T : Scalar + Serialize{

    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where S: Serializer {

        let mut state = serializer.serialize_struct("Cartessian", 2)?;
        state.serialize_field("dim", &N)?;
        state.serialize_field("coord", &Sequence(self.into_iter()))?;
        state.end()
    }
}


// ========================================================

// impl<T> Serialize for CartessianND<T>
//     where T : Scalar + Serialize{

//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//         where S: Serializer {

//         let mut state = serializer.serialize_struct("CartessianND", 2)?;
//         state.serialize_field("dim", &self.dim())?;
//         state.serialize_field("coord", &Sequence(self.into_iter()))?;
//         state.end()
//     }
// }



#[cfg(test)]
mod test {
    use crate::vector::CartessianND;

    #[test]
    fn test_serde_cartessian(){
        let a = CartessianND::new(vec![0, 3]);
        println!("{:?}", serde_json::to_string(&a).unwrap());
        println!("Hello");
    }
}
