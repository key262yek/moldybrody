use crate::boundary::plane::Parallelogram;
use crate::boundary::plane::PlanePair;
use crate::boundary::plane::SimpleBox;
use crate::boundary::plane::SimplePlanePair;
use crate::vector::Cartessian;
use crate::vector::Scalar;
use approx::AbsDiffEq;
use serde::de::{self, MapAccess, SeqAccess, Visitor};
use serde::{ser::SerializeStruct, Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::marker::PhantomData;

impl<T, const N: usize> Serialize for SimpleBox<T, N>
where
    [SimplePlanePair<T>; N]: Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("SimpleBox", 1)?;
        state.serialize_field("planes", &self.planes)?;
        state.end()
    }
}

struct SimpleBoxVisitor<T, const N: usize> {
    _marker_a: PhantomData<T>,
}

enum SimpleBoxField {
    Planes,
}

impl<T, const N: usize> SimpleBoxVisitor<T, N> {
    pub fn new() -> Self {
        SimpleBoxVisitor {
            _marker_a: PhantomData,
        }
    }
}

static SIMPLEBOX_FIELDS: &[&str] = &["planes"];

/// **Requires crate feature `"serde"`**
impl<'de, T, const N: usize> Deserialize<'de> for SimpleBox<T, N>
where
    [SimplePlanePair<T>; N]: Deserialize<'de>,
    T: Clone,
{
    fn deserialize<D>(deserializer: D) -> Result<SimpleBox<T, N>, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_struct("SimpleBox", SIMPLEBOX_FIELDS, SimpleBoxVisitor::new())
    }
}

impl<'de> Deserialize<'de> for SimpleBoxField {
    fn deserialize<D>(deserializer: D) -> Result<SimpleBoxField, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct SimpleBoxFieldVisitor;

        impl<'de> Visitor<'de> for SimpleBoxFieldVisitor {
            type Value = SimpleBoxField;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(r#""planes""#)
            }

            fn visit_str<E>(self, value: &str) -> Result<SimpleBoxField, E>
            where
                E: de::Error,
            {
                match value {
                    "planes" => Ok(SimpleBoxField::Planes),
                    other => Err(de::Error::unknown_field(other, SIMPLEBOX_FIELDS)),
                }
            }

            fn visit_bytes<E>(self, value: &[u8]) -> Result<SimpleBoxField, E>
            where
                E: de::Error,
            {
                match value {
                    b"planes" => Ok(SimpleBoxField::Planes),
                    other => Err(de::Error::unknown_field(
                        &format!("{:?}", other),
                        SIMPLEBOX_FIELDS,
                    )),
                }
            }
        }

        deserializer.deserialize_identifier(SimpleBoxFieldVisitor)
    }
}

impl<'de, T, const N: usize> Visitor<'de> for SimpleBoxVisitor<T, N>
where
    [SimplePlanePair<T>; N]: Deserialize<'de>,
    T: Clone,
{
    type Value = SimpleBox<T, N>;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("SimpleBox<T, N> representation")
    }

    fn visit_seq<V>(self, mut visitor: V) -> Result<SimpleBox<T, N>, V::Error>
    where
        V: SeqAccess<'de>,
    {
        let planes: [SimplePlanePair<T>; N] = match visitor.next_element()? {
            Some(value) => value,
            None => {
                return Err(de::Error::invalid_length(0, &self));
            }
        };

        Ok(SimpleBox::<T, N>::new(planes).unwrap())
    }

    fn visit_map<V>(self, mut visitor: V) -> Result<SimpleBox<T, N>, V::Error>
    where
        V: MapAccess<'de>,
    {
        let mut planes: Option<[SimplePlanePair<T>; N]> = None;

        while let Some(key) = visitor.next_key()? {
            match key {
                SimpleBoxField::Planes => {
                    planes = Some(visitor.next_value()?);
                }
            }
        }

        let planes = match planes {
            Some(planes) => planes,
            None => return Err(de::Error::missing_field("planes")),
        };

        Ok(SimpleBox::<T, N>::new(planes).unwrap())
    }
}

// ==================================================================================

impl<T, const N: usize> Serialize for Parallelogram<T, N>
where
    [PlanePair<Cartessian<T, N>>; N]: Serialize,
    T: Scalar + AbsDiffEq,
    <T as AbsDiffEq>::Epsilon: Clone,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Parallelogram", 1)?;
        state.serialize_field("planes", &self.planes)?;
        state.end()
    }
}

struct ParallelogramVisitor<T, const N: usize> {
    _marker_a: PhantomData<T>,
}

enum ParallelogramField {
    Planes,
}

impl<T, const N: usize> ParallelogramVisitor<T, N> {
    pub fn new() -> Self {
        ParallelogramVisitor {
            _marker_a: PhantomData,
        }
    }
}

static PARALLEOGRAM_FIELDS: &[&str] = &["planes"];

/// **Requires crate feature `"serde"`**
impl<'de, T, const N: usize> Deserialize<'de> for Parallelogram<T, N>
where
    [PlanePair<Cartessian<T, N>>; N]: Deserialize<'de>,
    T: Clone + Scalar + AbsDiffEq,
    <T as AbsDiffEq>::Epsilon: Clone,
{
    fn deserialize<D>(deserializer: D) -> Result<Parallelogram<T, N>, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_struct(
            "Parallelogram",
            PARALLEOGRAM_FIELDS,
            ParallelogramVisitor::new(),
        )
    }
}

impl<'de> Deserialize<'de> for ParallelogramField {
    fn deserialize<D>(deserializer: D) -> Result<ParallelogramField, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ParallelogramFieldVisitor;

        impl<'de> Visitor<'de> for ParallelogramFieldVisitor {
            type Value = ParallelogramField;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(r#""planes""#)
            }

            fn visit_str<E>(self, value: &str) -> Result<ParallelogramField, E>
            where
                E: de::Error,
            {
                match value {
                    "planes" => Ok(ParallelogramField::Planes),
                    other => Err(de::Error::unknown_field(other, PARALLEOGRAM_FIELDS)),
                }
            }

            fn visit_bytes<E>(self, value: &[u8]) -> Result<ParallelogramField, E>
            where
                E: de::Error,
            {
                match value {
                    b"planes" => Ok(ParallelogramField::Planes),
                    other => Err(de::Error::unknown_field(
                        &format!("{:?}", other),
                        PARALLEOGRAM_FIELDS,
                    )),
                }
            }
        }

        deserializer.deserialize_identifier(ParallelogramFieldVisitor)
    }
}

impl<'de, T, const N: usize> Visitor<'de> for ParallelogramVisitor<T, N>
where
    [PlanePair<Cartessian<T, N>>; N]: Deserialize<'de>,
    T: Clone + Scalar + AbsDiffEq,
    <T as AbsDiffEq>::Epsilon: Clone,
{
    type Value = Parallelogram<T, N>;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Parallelogram<T, N> representation")
    }

    fn visit_seq<V>(self, mut visitor: V) -> Result<Parallelogram<T, N>, V::Error>
    where
        V: SeqAccess<'de>,
    {
        let planes: [PlanePair<Cartessian<T, N>>; N] = match visitor.next_element()? {
            Some(value) => value,
            None => {
                return Err(de::Error::invalid_length(0, &self));
            }
        };

        Ok(Parallelogram::<T, N>::new(planes).unwrap())
    }

    fn visit_map<V>(self, mut visitor: V) -> Result<Parallelogram<T, N>, V::Error>
    where
        V: MapAccess<'de>,
    {
        let mut planes: Option<[PlanePair<Cartessian<T, N>>; N]> = None;

        while let Some(key) = visitor.next_key()? {
            match key {
                ParallelogramField::Planes => {
                    planes = Some(visitor.next_value()?);
                }
            }
        }

        let planes = match planes {
            Some(planes) => planes,
            None => return Err(de::Error::missing_field("planes")),
        };

        Ok(Parallelogram::<T, N>::new(planes).unwrap())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::vector::Cartessian;
    use serde_json::{from_str, to_string};

    #[test]
    fn test_serde_simple_box() {
        let simple_box: SimpleBox<f64, 2> = SimpleBox::from_pairs([[0.0, 1.0], [0.0, 1.0]]);
        let expected = r#"{"planes":[{"idx":0,"pos":[0.0,1.0]},{"idx":1,"pos":[0.0,1.0]}]}"#;
        assert_eq!(expected, to_string(&simple_box).unwrap());

        let expected: SimpleBox<f64, 2> = from_str(&expected).unwrap();
        assert_eq!(simple_box, expected);
    }

    #[test]
    fn test_serde_parallelogram() {
        let normal_vec1 = Cartessian::<f64, 2>::new([1.0, 0.0]);
        let const1 = [0.0f64, 1.0f64];
        let plane1: PlanePair<Cartessian<f64, 2>> = PlanePair::new(normal_vec1, const1).unwrap();

        let normal_vec2 = Cartessian::<f64, 2>::new([0.0, 1.0]);
        let const2 = [0.0f64, 1.0f64];
        let plane2: PlanePair<Cartessian<f64, 2>> = PlanePair::new(normal_vec2, const2).unwrap();

        let parallel: Parallelogram<f64, 2> = Parallelogram::new([plane1, plane2]).unwrap();
        let expected = r#"{"planes":[{"normal_vec":{"coord":[1.0,0.0]},"constant":[0.0,1.0]},{"normal_vec":{"coord":[0.0,1.0]},"constant":[0.0,1.0]}]}"#;
        assert_eq!(expected, to_string(&parallel).unwrap());

        let expected: Parallelogram<f64, 2> = from_str(&expected).unwrap();
        assert_eq!(parallel, expected);
    }
}
