//! Implementation of structures with the Topology trait corresponding specific point, neighbor pair
//!
//! [`Point`](trait@Point)와 [`Neighbor`](trait@Neighbor)에서 정의된 point, neighbor pair를 topology로 가지는 여러 시스템들을 정의하였습니다.
//! 점의 주변을 정의함으로써 점이 이동할 수 있는 영역을 확인한다던가, 가능한 이동인지 여부를 체크한다던가 하는 기능을 합니다.
//! 세부적인 기능은 [`Topology`](trait@Topology)에 설명되어 있습니다.

use crate::prelude::*;

/// Normal openball topology for continuous system
///
/// n차원 실공간의 topology는 보통 openball을 basis로 하여 주어집니다.
/// Continuous module은 연속적인 공간의 일반적 topology를 구현해 두었습니다.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct ContinuousTopology{
    /// Maximum step size available in the system
    max_step : f64,
}

impl Default for ContinuousTopology{
    fn default() -> Self{
        Self{
            max_step : 0f64,
        }
    }
}

impl ContinuousTopology{
    /// Initializing [`ContinuousTopology`](struct@ContinuousTopology) with a maximum step size.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// let sys = ContinuousTopology::new(0.1);
    /// ```
    pub fn new(max_step : f64) -> Self{
        Self{
            max_step,
        }
    }
}

impl<'a> Topology<'a, Cartessian1D<f64>, OpenBall1D<'a>> for ContinuousTopology{
    /// Return a openball of center pos.
    ///
    /// Cartessian1D 구조체로 주어지는 point를 중심으로 하는 OpenBall1D를 반환한다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// # use moldybrody::system::neighbor::OpenBall1D;
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let v = Cartessian1D{coord : 0.0};
    /// assert_eq!(sys.neighbor(&v), OpenBall1D::new(&v));
    /// ```
    fn neighbor<'b : 'a>(&self, pos : &'b Cartessian1D<f64>) -> OpenBall1D<'b>{
        OpenBall1D{
            center : &pos,
        }
    }

    /// Check whether the point is in a neighbor or not.
    ///
    /// OpenBall로 주어진 neighbor에 point가 속해있는지 여부를 확인해주는 함수이다.
    /// 시스템에는 OpenBall의 반지름이 max_step으로 정해져있는데,
    /// 이는 시뮬레이션 과정에서 입자가 움직일 수 있는 최대 변위를 계산하기 위해서다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// # use moldybrody::system::neighbor::OpenBall1D;
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let vt = Cartessian1D{coord : 0.05};
    ///
    /// let v1 = Cartessian1D{coord : 0.0};
    /// let n1 = OpenBall1D::new(&v1);
    /// assert_eq!(sys.inclusion(&vt, &n1), true);
    ///
    /// let v2 = Cartessian1D{coord : 1.0};
    /// let n2 = OpenBall1D::new(&v2);
    /// assert_eq!(sys.inclusion(&vt, &n2), false);
    /// ```
    ///
    /// # Panic
    ///
    /// point와 neighbor가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
    /// ```should_panic
    /// # use moldybrody::system::point::{Cartessian1D, Cartessian2D, CartessianND};
    /// # use moldybrody::system::neighbor::{OpenBall2D, OpenBallND};
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let v1 = Cartessian1D{coord : 0.0};
    ///
    /// let v2 = Cartessian2D{coord : [0.0; 2]};
    /// let n2 = OpenBall2D::new(&v2);
    /// // sys.inclusion(&v1, &n2);    // Cannot compile
    ///
    /// let v3 = CartessianND{coord : vec![0.0; 2]};
    /// let n3 = OpenBallND::new(&v3);
    /// // sys.inclusion(&v2, &n3);    // Cannot compile
    ///
    /// let v4 = CartessianND{coord : vec![0.0; 3]};
    /// let n4 = OpenBallND::new(&v4);
    /// sys.inclusion(&v3, &n4);       // Panic!
    /// ```
    fn inclusion(&self, pos : &Cartessian1D<f64>, neigh : &OpenBall1D) -> bool{
        let r = pos.distance(neigh.center);
        r < self.max_step
    }
}

impl<'a> Topology<'a, CartessianND<f64>, OpenBallND<'a>> for ContinuousTopology{
    /// Return a openball of center pos.
    ///
    /// CartessianND 구조체로 주어지는 point를 중심으로 하는 OpenBall1D를 반환한다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// # use moldybrody::system::neighbor::OpenBallND;
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let v = CartessianND{coord : vec![0.0; 3]};
    /// assert_eq!(sys.neighbor(&v), OpenBallND::new(&v));
    /// ```
    fn neighbor<'b : 'a>(&self, pos : &'b CartessianND<f64>) -> OpenBallND<'b>{
        OpenBallND{
            center : &pos,
        }
    }

    /// Check whether the point is in a neighbor or not.
    ///
    /// OpenBall로 주어진 neighbor에 point가 속해있는지 여부를 확인해주는 함수이다.
    /// 시스템에는 OpenBall의 반지름이 max_step으로 정해져있는데,
    /// 이는 시뮬레이션 과정에서 입자가 움직일 수 있는 최대 변위를 계산하기 위해서다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// # use moldybrody::system::neighbor::OpenBallND;
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let mut vt = CartessianND{coord : vec![0.0; 3]};
    /// vt.coord[0] = 0.05;
    ///
    /// let v1 = CartessianND{coord : vec![0.0; 3]};
    /// let n1 = OpenBallND::new(&v1);
    /// assert_eq!(sys.inclusion(&vt, &n1), true);
    ///
    /// let v2 = CartessianND{coord : vec![1.0; 3]};
    /// let n2 = OpenBallND::new(&v2);
    /// assert_eq!(sys.inclusion(&vt, &n2), false);
    /// ```
    ///
    /// # Panic
    ///
    /// point와 neighbor가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
    /// ```should_panic
    /// # use moldybrody::system::point::{Cartessian1D, Cartessian2D, CartessianND};
    /// # use moldybrody::system::neighbor::{OpenBall2D, OpenBallND};
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let v1 = Cartessian1D{coord : 0.0};
    ///
    /// let v2 = Cartessian2D{coord : [0.0; 2]};
    /// let n2 = OpenBall2D::new(&v2);
    /// // sys.inclusion(&v1, &n2);    // Cannot compile
    ///
    /// let v3 = CartessianND{coord : vec![0.0; 2]};
    /// let n3 = OpenBallND::new(&v3);
    /// // sys.inclusion(&v2, &n3);    // Cannot compile
    ///
    /// let v4 = CartessianND{coord : vec![0.0; 3]};
    /// let n4 = OpenBallND::new(&v4);
    /// sys.inclusion(&v3, &n4);       // Panic!
    /// ```
    fn inclusion(&self, pos : &CartessianND<f64>, neigh : &OpenBallND) -> bool{
        let r = pos.distance(neigh.center);
        r < self.max_step
    }
}


#[allow(unused_macros)]
macro_rules! impl_continuous_topology_nD{
    ($cartessian_name:ident, $neighbor_name:ident, $dim:expr) =>{
        impl<'a> Topology<'a, $cartessian_name<f64>, $neighbor_name<'a>> for ContinuousTopology{
            doc_comment!{
                concat!(
                    "Return a openball of center pos.

", stringify!($cartessian_name), " 구조체로 주어지는 point를 중심으로 하는 ", stringify!($neighbor_name), " 을 반환한다.

# Examples

```
# use moldybrody::system::point::", stringify!($cartessian_name), ";
# use moldybrody::system::neighbor::", stringify!($neighbor_name), ";
# use moldybrody::system::topology::ContinuousTopology;
# use moldybrody::system::Topology;let sys = ContinuousTopology::new(0.1);
let v = ", stringify!($cartessian_name), "{coord : [0.0; ", $dim, "]};
assert_eq!(sys.neighbor(&v), ", stringify!($neighbor_name), "::new(&v));
```"
                ),
                fn neighbor<'b : 'a>(&self, pos : &'b $cartessian_name<f64>) -> $neighbor_name<'b>{
                    $neighbor_name{
                        center : &pos,
                    }
                }
            }

            doc_comment!{
                concat!(
                    "Check whether the point is in a neighbor or not.

OpenBall로 주어진 neighbor에 point가 속해있는지 여부를 확인해주는 함수이다.
시스템에는 OpenBall의 반지름이 max_step으로 정해져있는데,
이는 시뮬레이션 과정에서 입자가 움직일 수 있는 최대 변위를 계산하기 위해서다.

# Examples

```
# use moldybrody::system::point::", stringify!($cartessian_name), ";
# use moldybrody::system::neighbor::", stringify!($neighbor_name), ";
# use moldybrody::system::topology::ContinuousTopology;
# use moldybrody::system::Topology;
let sys = ContinuousTopology::new(0.1);
let mut vt = ", stringify!($cartessian_name), "{coord : [0.0; ", $dim, "]};
vt.coord[0] = 0.05;

let v1 = ", stringify!($cartessian_name), "{coord : [0.0; ", $dim, "]};
let n1 = ", stringify!($neighbor_name), "::new(&v1);
assert_eq!(sys.inclusion(&vt, &n1), true);

let v2 = ", stringify!($cartessian_name), "{coord : [1.0; ", $dim, "]};
let n2 = ", stringify!($neighbor_name), "::new(&v2);
assert_eq!(sys.inclusion(&vt, &n2), false);
```

# Panic

point와 neighbor가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
```should_panic
# use moldybrody::system::point::{Cartessian1D, Cartessian2D, CartessianND};
# use moldybrody::system::neighbor::{OpenBall2D, OpenBallND};
# use moldybrody::system::topology::ContinuousTopology;
# use moldybrody::system::Topology;
let sys = ContinuousTopology::new(0.1);
let v1 = Cartessian1D{coord : 0.0};

let v2 = Cartessian2D{coord : [0.0; 2]};
let n2 = OpenBall2D::new(&v2);
// sys.inclusion(&v1, &n2);    // Cannot compile

let v3 = CartessianND{coord : vec![0.0; 2]};
let n3 = OpenBallND::new(&v3);
// sys.inclusion(&v2, &n3);    // Cannot compile

let v4 = CartessianND{coord : vec![0.0; 3]};
let n4 = OpenBallND::new(&v4);
sys.inclusion(&v3, &n4);       // Panic!
```"
                ),
                fn inclusion(&self, pos : &$cartessian_name<f64>, neigh : &$neighbor_name) -> bool{
                    let r = pos.distance(neigh.center);
                    r < self.max_step
                }
            }
        }
    }
}

impl_continuous_topology_nD!(Cartessian2D, OpenBall2D, 2);
impl_continuous_topology_nD!(Cartessian3D, OpenBall3D, 3);
impl_continuous_topology_nD!(Cartessian4D, OpenBall4D, 4);

