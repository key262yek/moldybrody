use super::Cartessian;
use serde::de::{self, MapAccess, SeqAccess, Visitor};
use serde::{ser::SerializeStruct, Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::marker::PhantomData;

impl<T, const N: usize> Serialize for Cartessian<T, N>
where
    [T; N]: Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Cartessian", 1)?;
        state.serialize_field("coord", &self.coord)?;
        state.end()
    }
}

struct CartessianVisitor<T, const N: usize> {
    _marker_a: PhantomData<T>,
}

enum CartessianField {
    Coord,
}

impl<T, const N: usize> CartessianVisitor<T, N> {
    pub fn new() -> Self {
        CartessianVisitor {
            _marker_a: PhantomData,
        }
    }
}

static CARTESSIAN_FIELDS: &[&str] = &["coord"];

/// **Requires crate feature `"serde"`**
impl<'de, T, const N: usize> Deserialize<'de> for Cartessian<T, N>
where
    [T; N]: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Cartessian<T, N>, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_struct("Cartessian", CARTESSIAN_FIELDS, CartessianVisitor::new())
    }
}

impl<'de> Deserialize<'de> for CartessianField {
    fn deserialize<D>(deserializer: D) -> Result<CartessianField, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct CartessianFieldVisitor;

        impl<'de> Visitor<'de> for CartessianFieldVisitor {
            type Value = CartessianField;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(r#""coord""#)
            }

            fn visit_str<E>(self, value: &str) -> Result<CartessianField, E>
            where
                E: de::Error,
            {
                match value {
                    "coord" => Ok(CartessianField::Coord),
                    other => Err(de::Error::unknown_field(other, CARTESSIAN_FIELDS)),
                }
            }

            fn visit_bytes<E>(self, value: &[u8]) -> Result<CartessianField, E>
            where
                E: de::Error,
            {
                match value {
                    b"coord" => Ok(CartessianField::Coord),
                    other => Err(de::Error::unknown_field(
                        &format!("{:?}", other),
                        CARTESSIAN_FIELDS,
                    )),
                }
            }
        }

        deserializer.deserialize_identifier(CartessianFieldVisitor)
    }
}

impl<'de, T, const N: usize> Visitor<'de> for CartessianVisitor<T, N>
where
    [T; N]: Deserialize<'de>,
{
    type Value = Cartessian<T, N>;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Cartessian<T, N> representation")
    }

    fn visit_seq<V>(self, mut visitor: V) -> Result<Cartessian<T, N>, V::Error>
    where
        V: SeqAccess<'de>,
    {
        let coord: [T; N] = match visitor.next_element()? {
            Some(value) => value,
            None => {
                return Err(de::Error::invalid_length(0, &self));
            }
        };

        Ok(Cartessian::<T, N>::new(coord))
    }

    fn visit_map<V>(self, mut visitor: V) -> Result<Cartessian<T, N>, V::Error>
    where
        V: MapAccess<'de>,
    {
        let mut coord: Option<[T; N]> = None;

        while let Some(key) = visitor.next_key()? {
            match key {
                CartessianField::Coord => {
                    coord = Some(visitor.next_value()?);
                }
            }
        }

        let coord = match coord {
            Some(coord) => coord,
            None => return Err(de::Error::missing_field("coord")),
        };

        Ok(Cartessian::<T, N>::new(coord))
    }
}

#[cfg(test)]
mod test {
    use crate::vector::{Cartessian, CartessianND};
    use serde_json::{from_str, to_string};

    #[test]
    fn test_serde_cartessian() {
        let a: Cartessian<i32, 2> = Cartessian::new([0, 3]);
        let expected = r#"{"coord":[0,3]}"#;
        assert_eq!(expected, to_string(&a).unwrap());

        let expected: Cartessian<i32, 2> = from_str(&expected).unwrap();
        assert_eq!(a, expected);
    }

    #[test]
    fn test_serde_cartessian_nd() {
        let a: CartessianND<i32> = CartessianND::new(vec![0, 3]);
        let expected = r#"{"coord":[0,3]}"#;
        assert_eq!(expected, to_string(&a).unwrap());

        let expected: CartessianND<i32> = from_str(&expected).unwrap();
        assert_eq!(a, expected);
    }
}
