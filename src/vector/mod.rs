use approx::AbsDiffEq;
use ndarray_linalg::types::Scalar as NDScalar;
use num_traits::{One, Zero};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};

use std::fmt::{Debug, Display};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

pub trait Integer:
    'static
    + Copy
    + Eq
    + PartialEq
    + PartialOrd
    + Ord
    + Zero
    + One
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + AddAssign
    + SubAssign
    + DivAssign
    + MulAssign
{
}
pub trait Scalar:
    'static
    + Debug
    + Copy
    + PartialEq
    + Sum
    + Zero
    + One
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + AddAssign
    + SubAssign
    + DivAssign
    + MulAssign
    + Serialize
    + DeserializeOwned
{
}
pub trait Float: NDScalar + Display + Debug + Scalar + Send + Sync + AbsDiffEq {}

pub mod arithmetic;
pub mod basic;
pub mod curvilinear;
pub mod product;
pub mod vector_serde;

pub trait Vector {
    type Item: Scalar;
    /// Return dimension of vector
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a : Cartessian<f64, 2> = Cartessian::new([3.0, 5.0]);
    /// assert_eq!(a.dim(), 2);
    /// ```
    fn dim(&self) -> usize;
}

impl<T: Scalar, const N: usize> Vector for [T; N] {
    type Item = T;

    fn dim(&self) -> usize {
        N
    }
}

pub trait Dim<const N: usize> {}
impl<T, const N: usize> Dim<N> for Cartessian<T, N> {}
impl<T> Dim<0> for CartessianND<T> {}
impl<T> Dim<1> for CartessianND<T> {}
impl<T> Dim<2> for CartessianND<T> {}
impl<T> Dim<3> for CartessianND<T> {}
impl<T> Dim<4> for CartessianND<T> {}
impl<T> Dim<5> for CartessianND<T> {}

impl<V: Vector, W: Vector> Vector for (V, W) {
    type Item = V::Item;

    fn dim(&self) -> usize {
        self.0.dim() + self.1.dim()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
/// Structrure wrap an array as coordinate
///
/// Array를 coordinate로 삼은 벡터 구조체.
/// Array로 만들어졌기 때문에 컴파일 과정에서 서로 차원이 다른 두 벡터가 연산될 가능성이 차단된다.
pub struct Cartessian<T, const N: usize> {
    pub coord: [T; N],
}

pub type Cartessian1D<T> = Cartessian<T, 1>;
pub type Cartessian2D<T> = Cartessian<T, 2>;
pub type Cartessian3D<T> = Cartessian<T, 3>;
pub type Cartessian4D<T> = Cartessian<T, 4>;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
/// Structrure wrap an vector as coordinate
///
/// Vector를 coordinate로 삼은 벡터 구조체.
/// 벡터 길이가 자유롭게 선언될 수 있단 측면에서 시뮬레이션 시스템의 차원을 변경하기 수월하지만,
/// 매 연산마다 벡터 길이가 다른지 여부를 확인해야한다는 단점이 있다.
pub struct CartessianND<T> {
    pub coord: Vec<T>,
}

impl<T, const N: usize> Cartessian<T, N> {
    /// Create a new vector
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a : Cartessian<f64, 2> = Cartessian::new([2.0, 3.0]);
    /// ```
    pub fn new(coord: [T; N]) -> Self {
        Self { coord: coord }
    }

    /// Create a new vector from Vec<T>
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a : Cartessian<f64, 2> = Cartessian::from_vec(vec![2.0, 3.0]).unwrap();
    /// ```
    // pub fn from_vec(coord: Vec<T>) -> Result<Self, Error> {
    //     if coord.len() != N {
    //         return Err(Error::make_error_syntax(ErrorCode::InvalidDimension));
    //     }

    //     Ok(Self {
    //         coord: coord.try_into().unwrap(),
    //     })
    // }

    /// Create a new vector from Iterator<Item = T>
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let iter = vec![2.0, 3.0].into_iter().map(|x| x * x);
    /// let a : Cartessian<f64, 2> = Cartessian::from_iter(iter).unwrap();
    /// ```
    // pub fn from_iter<I>(iterable: I) -> Result<Self, Error>
    // where
    //     I: IntoIterator<Item = T>,
    // {
    //     Self::from_vec(iterable.into_iter().collect::<Vec<T>>())
    // }

    // #[allow(dead_code)]
    // pub(crate) fn from_vec_unchecked(coord: Vec<T>) -> Self {
    //     Self {
    //         coord: coord.try_into().unwrap(),
    //     }
    // }

    // #[allow(dead_code)]
    // pub(crate) fn from_iter_unchecked<I>(iterable: I) -> Self
    // where
    //     I: IntoIterator<Item = T>,
    // {
    //     Self::from_vec_unchecked(iterable.into_iter().collect::<Vec<T>>())
    // }

    /// Return length of coordinate
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a : Cartessian<f64, 2> = Cartessian::new([3.0, 5.0]);
    /// assert_eq!(a.len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        N
    }

    pub fn dim(&self) -> usize {
        N
    }
}

impl<T, const N: usize> Vector for Cartessian<T, N>
where
    T: Scalar,
{
    type Item = T;
    /// Return dimension of vector
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a : Cartessian<f64, 2> = Cartessian::new([3.0, 5.0]);
    /// assert_eq!(a.dim(), 2);
    /// ```
    fn dim(&self) -> usize {
        N
    }
}

impl<T, const N: usize> Vector for &Cartessian<T, N>
where
    T: Scalar,
{
    type Item = T;
    /// Return dimension of vector
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a : Cartessian<f64, 2> = Cartessian::new([3.0, 5.0]);
    /// assert_eq!(a.dim(), 2);
    /// ```
    fn dim(&self) -> usize {
        N
    }
}

impl<T, const N: usize> Vector for &mut Cartessian<T, N>
where
    T: Scalar,
{
    type Item = T;
    /// Return dimension of vector
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a : Cartessian<f64, 2> = Cartessian::new([3.0, 5.0]);
    /// assert_eq!(a.dim(), 2);
    /// ```
    fn dim(&self) -> usize {
        N
    }
}

impl<T> CartessianND<T> {
    /// Create a new vector
    ///
    /// ```
    /// # use moldybrody::vector::CartessianND;
    /// let a : CartessianND<f64> = CartessianND::new(vec![2.0, 3.0]);
    /// ```
    pub fn new(coord: Vec<T>) -> Self {
        Self { coord: coord }
    }

    /// Create a new vector from Iterator<Item = T>
    ///
    /// ```
    /// # use moldybrody::vector::CartessianND;
    /// let iter = vec![2.0, 3.0].into_iter().map(|x| x * x);
    /// let a : CartessianND<f64> = CartessianND::from_iter(iter);
    /// ```
    pub fn from_iter<I>(iterable: I) -> Self
    where
        I: IntoIterator<Item = T>,
    {
        Self::new(iterable.into_iter().collect::<Vec<T>>())
    }

    /// Return length of coordinate
    ///
    /// ```
    /// # use moldybrody::vector::CartessianND;
    /// let a : CartessianND<f64> = CartessianND::new(vec![3.0, 5.0]);
    /// assert_eq!(a.len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.coord.len()
    }

    pub fn dim(&self) -> usize {
        self.coord.len()
    }
}

impl<T: Scalar> Vector for CartessianND<T> {
    type Item = T;
    /// Return dimension of vector
    ///
    /// ```
    /// # use moldybrody::vector::CartessianND;
    /// let a : CartessianND<f64> = CartessianND::new(vec![3.0, 5.0]);
    /// assert_eq!(a.dim(), 2);
    /// ```
    fn dim(&self) -> usize {
        self.coord.len()
    }
}

impl<T: Scalar> Vector for &CartessianND<T> {
    type Item = T;
    /// Return dimension of vector
    ///
    /// ```
    /// # use moldybrody::vector::CartessianND;
    /// let a : CartessianND<f64> = CartessianND::new(vec![3.0, 5.0]);
    /// assert_eq!(a.dim(), 2);
    /// ```
    fn dim(&self) -> usize {
        self.coord.len()
    }
}

impl<T: Scalar> Vector for &mut CartessianND<T> {
    type Item = T;
    /// Return dimension of vector
    ///
    /// ```
    /// # use moldybrody::vector::CartessianND;
    /// let a : CartessianND<f64> = CartessianND::new(vec![3.0, 5.0]);
    /// assert_eq!(a.dim(), 2);
    /// ```
    fn dim(&self) -> usize {
        self.coord.len()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_new() {
        let a = Cartessian::new([2; 3]);
        assert_eq!(a.coord, [2; 3]);

        // let a = Cartessian::<i32, 3>::from_vec(vec![2; 3]).unwrap();
        // assert_eq!(a.coord, [2; 3]);

        // let a = Cartessian::<i32, 3>::from_iter(vec![2; 3].into_iter()).unwrap();
        // assert_eq!(a.coord, [2; 3]);

        assert_eq!(a.len(), 3);
        assert_eq!(a.dim(), 3);
    }

    #[test]
    fn test_new_nd() {
        let a = CartessianND::new(vec![1f64, 0f64]);
        assert_eq!(a.coord, vec![1f64, 0f64]);

        let a = CartessianND::from_iter(vec![2; 3].into_iter());
        assert_eq!(a.coord, [2; 3]);

        assert_eq!(a.len(), 3);
        assert_eq!(a.dim(), 3);
    }
}
