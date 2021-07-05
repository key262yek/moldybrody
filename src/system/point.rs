//! Specific implementation for structures with the Point trait
//!
//! 여러 시스템에서 Point 역할을 할 구조체들을 정의하고, 기능들을 포함하고 있습니다.

use crate::prelude::*;

/// Cartessian coordinate of 1-dimensional continuous system
///
/// 1D의 경우는 굳이 Iteraion을 하지 않아도 되는 경우입니다.
/// 비효율을 줄이기 위해 따로 정의합니다.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct Cartessian1D<T>{
    /// Coordinate of Point
    pub coord : T,
}

impl Point for Cartessian1D<f64>{}

impl Cartessian1D<f64>{
    /// return a dimension of system
    ///
    /// 시스템의 차원은 coordinate vector의 length로 주어진다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// let v = Cartessian1D{coord : 0.0};
    /// assert_eq!(v.dim(), 1);
    /// ```
    pub fn dim(&self) -> usize{
        1
    }

    /// return a distance to other point
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// let v1 = Cartessian1D{coord : 0.0};
    /// let mut v2 = Cartessian1D{coord : 1.0};
    ///
    /// assert_eq!(v1.distance(&v2), 1.0);
    /// ```
    pub fn distance(&self, other : &Self) -> f64{
        (self.coord - other.coord).abs()
    }

    /// return a norm of vector
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// let v1 = Cartessian1D{coord : 0.0};
    /// assert_eq!(v1.norm(), 0.0);
    ///
    /// let v2 = Cartessian1D{coord : 3.0};
    /// assert_eq!(v2.norm(), 3.0);
    ///
    /// let v3 = Cartessian1D{coord : -3.0};
    /// assert_eq!(v3.norm(), 3.0);
    /// ```
    pub fn norm(&self) -> f64{
        return self.coord.abs();
    }
}

impl Point for Cartessian1D<i32>{}

impl Cartessian1D<i32>{
    /// return a dimension of point
    ///
    /// 시스템의 차원은 coordinate vector의 length로 주어진다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// let v = Cartessian1D{coord : 0};
    /// assert_eq!(v.dim(), 1);
    /// ```
    pub fn dim(&self) -> usize{
        1
    }

    /// return a taxi distance to other point
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// let v1 = Cartessian1D{coord : 0};
    /// let mut v2 = Cartessian1D{coord : 1};
    ///
    /// assert_eq!(v1.taxi_distance(&v2), 1);
    /// ```
    pub fn taxi_distance(&self, other : &Self) -> usize{
        (self.coord - other.coord).abs() as usize
    }

    /// return a taxi norm of vector
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// let v1 = Cartessian1D{coord : 0};
    /// assert_eq!(v1.taxi_norm(), 0);
    ///
    /// let v2 = Cartessian1D{coord : 3};
    /// assert_eq!(v2.taxi_norm(), 3);
    ///
    /// let v3 = Cartessian1D{coord : -3};
    /// assert_eq!(v2.taxi_norm(), 3);
    /// ```
    pub fn taxi_norm(&self) -> usize{
        return self.coord.abs() as usize;
    }
}




/// Cartessian coordinate of n-dimensional continuous system
///
/// [`Cartessian1D`](struct@Cartessian1D) 같이 fixed sized array를 coordinate로 하면 비효율은 줄어들 수 있으나,
/// 일반적인 차원의 시스템을 모두 정의하는 것은 불가능해집니다.
/// 이에 fixed sized array가 아닌 vector를 이용해 일반 차원의 좌표계를 정의했습니다.
#[derive(Clone, Debug, PartialEq, PartialOrd, Hash)]
pub struct CartessianND<T>{
    /// Coordinate of Point
    pub coord : Vec<T>,
}

impl Point for CartessianND<f64>{}

impl CartessianND<f64>{
    /// return a dimension of system
    ///
    /// 시스템의 차원은 coordinate vector의 length로 주어진다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// let v = CartessianND{coord : vec![0.0; 3]};
    /// assert_eq!(v.dim(), 3);
    /// ```
    pub fn dim(&self) -> usize{
        self.coord.len()
    }

    /// return a distance to other point
    /// 
    /// # Examples 
    /// 
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// let v1 = CartessianND{coord : vec![0.0; 2]};
    /// let v2 = CartessianND{coord : vec![3.0, 4.0]};
    /// assert_eq!(v1.distance(&v2), 5.0);
    /// ```
    ///
    /// # Panic
    ///
    /// 이 함수는 두 point 사이 거리를 반환하는데, 만약 두 point의 dimension이 다르면 panic한다.
    /// ```should_panic
    /// # use moldybrody::system::point::CartessianND;
    /// let v1 = CartessianND{coord : vec![0.0; 3]};
    /// let v2 = CartessianND{coord : vec![3.0, 4.0]};
    /// let r = v1.distance(&v2);
    /// ```
    pub fn distance(&self, other : &Self) -> f64{
        let dim = self.dim();
        if dim != other.dim(){
            panic!("{}", ErrorCode::InvalidDimension);
        }

        let mut r = 0f64;
        for i in 0..dim{
            let dx = self.coord[i] - other.coord[i];
            r += dx * dx;
        }

        return r.sqrt();
    }

    /// return a norm of vector
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// let v1 = CartessianND{coord : vec![0.0; 2]};
    /// assert_eq!(v1.norm(), 0.0);
    ///
    /// let v2 = CartessianND{coord : vec![3.0, 4.0]};
    /// assert_eq!(v2.norm(), 5.0);
    /// ```
    pub fn norm(&self) -> f64{
        let mut r = 0f64;

        for x in &self.coord{
            r += x * x;
        }

        return r.sqrt();
    }
}

impl Point for CartessianND<i32>{}

impl CartessianND<i32>{
    /// return a dimension of point
    ///
    /// 시스템의 차원은 coordinate vector의 length로 주어진다.
    /// 
    /// # Examples 
    /// 
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// let v = CartessianND{coord : vec![0; 3]};
    /// assert_eq!(v.dim(), 3);
    /// ```
    pub fn dim(&self) -> usize{
        self.coord.len()
    }

    /// return a taxi distance to other point
    /// 
    /// # Examples 
    /// 
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// let v1 = CartessianND{coord : vec![0; 2]};
    /// let v2 = CartessianND{coord : vec![3, 4]};
    /// assert_eq!(v1.taxi_distance(&v2), 7);
    /// ```
    ///
    /// # Panic
    ///
    /// 이 함수는 두 point 사이 거리를 반환하는데, 만약 두 point의 dimension이 다르면 panic한다.
    /// ```should_panic
    /// # use moldybrody::system::point::CartessianND;
    /// let v1 = CartessianND{coord : vec![0; 3]};
    /// let v2 = CartessianND{coord : vec![3, 4]};
    /// let r = v1.taxi_distance(&v2);
    /// ```
    pub fn taxi_distance(&self, other : &Self) -> usize{
        let dim = self.dim();
        if dim != other.dim(){
            panic!("{}", ErrorCode::InvalidDimension);
        }

        let mut r = 0usize;
        for i in 0..dim{
            let dx = self.coord[i] - other.coord[i] ;
            r += dx.abs() as usize;
        }
        return r;
    }

    /// return a taxi norm of vector
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// let v1 = CartessianND{coord : vec![0; 2]};
    /// assert_eq!(v1.taxi_norm(), 0);
    ///
    /// let v2 = CartessianND{coord : vec![3, 4]};
    /// assert_eq!(v2.taxi_norm(), 7);
    ///
    /// let v3 = CartessianND{coord : vec![3, -4]};
    /// assert_eq!(v3.taxi_norm(), 7);
    /// ```
    pub fn taxi_norm(&self) -> usize{
        let mut r = 0;
        for x in &self.coord{
            r += x.abs()
        }
        return r as usize;
    }
}


#[allow(unused_macros)]
macro_rules! impl_cartessian_nD{
    ($type_name:ident, $dim:expr) =>{
        doc_comment!{
           concat!(
                "Cartessian coordinate of ", $dim, "-dimensional continuous system.",
"

Fixed sized array를 coordinate로 하는 cartessian 좌표계입니다.
이렇게 구조체를 선언하면 서로 다른 dimension의 좌표계와 같이 계산될 가능성이 원천 차단되어 compile 단계에서 많은 오류를 배제할 수 있고,
또한 계산 전 dimension을 체크하는 비효율을 줄일 수 있습니다. 게다가 size가 정해져있어 Copy가 가능하단 점 역시 매우 유용한 점입니다."
           ),
            #[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
            pub struct $type_name<T>{
                /// Coordinate of Point
                pub coord : [T; $dim],
            }
        }

        impl Point for $type_name<f64>{}

        impl $type_name<f64>{
            doc_comment!{
                concat!(
                    "return a dimension of system.

시스템의 차원은 coordinate vector의 length로 주어진다.

# Examples

```
# use moldybrody::system::point::", stringify!($type_name), ";
let v = ", stringify!($type_name), "{coord : [0.0; ", $dim, "]};
assert_eq!(v.dim(), ", $dim, ");
```"
                ),
                pub fn dim(&self) -> usize{
                    $dim
                }
            }

            doc_comment!{
                concat!(
                    "return a distance to other point

# Examples

```
# use moldybrody::system::point::", stringify!($type_name), ";
let v1 = ", stringify!($type_name), "{coord : [0.0; ", $dim, "]};
let mut v2 = ", stringify!($type_name), "{coord : [0.0; ", $dim, "]};
v2.coord[0] = 1.0;

assert_eq!(v1.distance(&v2), 1.0);
```"),
                pub fn distance(&self, other : &Self) -> f64{
                    let mut r = 0f64;
                    for i in 0..$dim{
                        let dx = self.coord[i] - other.coord[i];
                        r += dx * dx;
                    }

                    return r.sqrt();
                }
            }

            doc_comment!{
                concat!(
                    "return a norm of vector

# Examples

```
# use moldybrody::system::point::", stringify!($type_name), ";
let mut v1 = ", stringify!($type_name), "{coord : [0.0; ", $dim, "]};
assert_eq!(v1.norm(), 0.0);

v1.coord[0] = 3.0;
assert_eq!(v1.norm(), 3.0);

v1.coord[1] = -4.0;
assert_eq!(v1.norm(), 5.0);
```"),
                pub fn norm(&self) -> f64{
                    let mut r = 0f64;
                    for x in &self.coord{
                        r += x * x;
                    }

                    return r.sqrt();
                }
            }
        }

        impl Point for $type_name<i32>{}

        impl $type_name<i32>{
            doc_comment!{
                concat!(
                    "return a dimension of point.

시스템의 차원은 coordinate vector의 length로 주어진다.

# Examples

```
# use moldybrody::system::point::", stringify!($type_name), ";
let v = ", stringify!($type_name), "{coord : [0; ", $dim, "]};
assert_eq!(v.dim(), ", $dim, ");
```"),
                pub fn dim(&self) -> usize{
                    $dim
                }
            }

            doc_comment!{
                concat!(
                    "return a taxi distance to other point.

# Examples

```
# use moldybrody::system::point::", stringify!($type_name),";
let v1 = ", stringify!($type_name),"{coord : [0; ", $dim, "]};
let mut v2 = ", stringify!($type_name),"{coord : [0; ", $dim, "]};
v2.coord[0] = 1;
///
assert_eq!(v1.taxi_distance(&v2), 1);
```
"),
                pub fn taxi_distance(&self, other : &Self) -> usize{
                    let mut r = 0usize;
                    for i in 0..$dim{
                        let dx = self.coord[i] - other.coord[i] ;
                        r += dx.abs() as usize;
                    }
                    return r;
                }
            }

            doc_comment!{
                concat!(
                    "return a taxi norm of vector

# Examples

```
# use moldybrody::system::point::", stringify!($type_name), ";
let mut v1 = ", stringify!($type_name), "{coord : [0; ", $dim, "]};
assert_eq!(v1.taxi_norm(), 0);

v1.coord[0] = 3;
assert_eq!(v1.taxi_norm(), 3);

v1.coord[1] = -4;
assert_eq!(v1.taxi_norm(), 7);
```"
                ),
                pub fn taxi_norm(&self) -> usize{
                    let mut r = 0;
                    for x in &self.coord{
                        r += x.abs()
                    }
                    return r as usize;
                }
            }
        }
    }
}

impl_cartessian_nD!(Cartessian2D, 2);
impl_cartessian_nD!(Cartessian3D, 3);
impl_cartessian_nD!(Cartessian4D, 4);


// ===========================================================================================
// ===========================================================================================
// ===========================================================================================


/// node of a network
///
/// Network에서 node는 index만 부여된 상태로 이해하면 됩니다.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct NodeIndex{
    /// Index of node
    pub index : usize,
}

impl Point for NodeIndex{}

