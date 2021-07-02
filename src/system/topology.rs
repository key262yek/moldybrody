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

impl Topology<Cartessian1D<f64>> for ContinuousTopology{
    /// Check whether the movement is valid or not
    ///
    /// point의 이동이 가능한 move인지 아닌지 여부를 확인해주는 함수이다.
    /// 시스템은 max_step이 정해져있는데,
    /// 이는 시뮬레이션 과정에서 입자가 움직일 수 있는 최대 변위를 제한한다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let vt = Cartessian1D{coord : 0.0};
    ///
    /// let move1 = Cartessian1D{coord : 0.05};
    /// assert_eq!(sys.check_move(&vt, &move1), true);
    ///
    /// let move2 = Cartessian1D{coord : 0.11};
    /// assert_eq!(sys.check_move(&vt, &move2), false);
    /// ```
    ///
    /// # Panic
    ///
    /// point와 movement가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
    /// ```should_panic
    /// # use moldybrody::system::point::{Cartessian1D, Cartessian2D, CartessianND};
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let p = Cartessian2D{coord : [1.0; 2]};
    /// let move = CartessianND{coord : vec![0.0; 2]};
    /// // sys.check_move(&p, &move);    // Cannot compile
    ///
    /// let p2 = CartessianND{coord : vec![0.0; 2]};
    /// let move2 = CartessianND{coord : vec![0.0; 3]};
    /// sys.check_move(&p2, &move2);       // Panic!
    /// ```
    fn check_move(&self, _pos : &Cartessian1D<f64>, movement : &Cartessian1D<f64>) -> bool{
        let r = movement.norm();
        r < self.max_step
    }
}

impl Topology<CartessianND<f64>> for ContinuousTopology{
    /// Check whether the movement is valid or not
    ///
    /// point의 이동이 가능한 move인지 아닌지 여부를 확인해주는 함수이다.
    /// 시스템은 max_step이 정해져있는데,
    /// 이는 시뮬레이션 과정에서 입자가 움직일 수 있는 최대 변위를 제한한다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let vt = CartessianND{coord : vec![2.0; 2]};
    ///
    /// let mut move1 = CartessianND{coord : vec![0.0, 0.05]};   // length of movement is smaller than 0.1
    /// assert_eq!(sys.check_move(&vt, &move1), true);
    ///
    /// let move2 = CartessianND{coord : vec![0.1, 0.05]};       // length of movement is larger than 0.1
    /// assert_eq!(sys.check_move(&vt, &move2), false);
    /// ```
    ///
    /// # Panic
    ///
    /// point와 movement가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
    /// ```should_panic
    /// # use moldybrody::system::point::{Cartessian1D, Cartessian2D, CartessianND};
    /// # use moldybrody::system::topology::ContinuousTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = ContinuousTopology::new(0.1);
    /// let p = Cartessian2D{coord : [1.0; 2]};
    /// let move = CartessianND{coord : vec![0.0; 2]};
    /// // sys.check_move(&p, &move);    // Cannot compile
    ///
    /// let p2 = CartessianND{coord : vec![0.0; 2]};
    /// let move2 = CartessianND{coord : vec![0.0; 3]};
    /// sys.check_move(&p2, &move2);       // Panic!
    /// ```
    fn check_move(&self, _pos : &CartessianND<f64>, movement : &CartessianND<f64>) -> bool{
        let r = movement.norm();
        r < self.max_step
    }
}


#[allow(unused_macros)]
macro_rules! impl_continuous_topology_nD{
    ($cartessian_name:ident, $dim:expr) =>{
        impl Topology<$cartessian_name<f64>> for ContinuousTopology{

            doc_comment!{
                concat!(
                    "Check whether the movement is valid or not

point의 이동이 가능한 move인지 아닌지 여부를 확인해주는 함수이다.
시스템은 max_step이 정해져있는데,
이는 시뮬레이션 과정에서 입자가 움직일 수 있는 최대 변위를 제한한다.

# Examples

```
# use moldybrody::system::point::", stringify!($cartessian_name), ";
# use moldybrody::system::topology::ContinuousTopology;
# use moldybrody::system::Topology;
let sys = ContinuousTopology::new(0.1);
let vt = ", stringify!($cartessian_name), "{coord : [2.0; ", $dim, "]};
let mut move1 = ", stringify!($cartessian_name), "{coord : [0.0; 2]};

move1.coord[0] = 0.05;   // length of movement is smaller than 0.1
assert_eq!(sys.check_move(&vt, &move1), true);

move1.coord[1] = 0.1;       // length of movement is larger than 0.1
assert_eq!(sys.check_move(&vt, &move1), false);
```

# Panic

point와 movement가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
```should_panic
# use moldybrody::system::point::{Cartessian1D, Cartessian2D, CartessianND};
# use moldybrody::system::topology::ContinuousTopology;
# use moldybrody::system::Topology;
let sys = ContinuousTopology::new(0.1);
let p = Cartessian2D{coord : [1.0; 2]};
let move = CartessianND{coord : vec![0.0; 2]};
// sys.check_move(&p, &move);    // Cannot compile

let p2 = CartessianND{coord : vec![0.0; 2]};
let move2 = CartessianND{coord : vec![0.0; 3]};
sys.check_move(&p2, &move2);       // Panic!
```"
                ),
                fn check_move(&self, _pos : &$cartessian_name<f64>, movement : &$cartessian_name<f64>) -> bool{
                    let r = movement.norm();
                    r < self.max_step
                }
            }
        }
    }
}

impl_continuous_topology_nD!(Cartessian2D, 2);
impl_continuous_topology_nD!(Cartessian3D, 3);
impl_continuous_topology_nD!(Cartessian4D, 4);

