//! Implementation of structures with the Topology trait corresponding specific point, neighbor pair
//!
//! [`Point`](trait@Point)에 정의된 여러 structure들의 topology를 정의하였습니다.
//! 점이 특정한 점으로 이동할 수 있는지 체크하는 기능을 합니다.
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

// ===========================================================================================
// ===========================================================================================
// ===========================================================================================

/// Normal near neighbor topology for discrete system
///
/// n차원 lattice의 topology는 lattice의 주변 node들로 주어집니다.
/// 한 번에 한 칸만 갈 수 있는 nearest neighbor 경우도 있지만,
/// 한 번에 두세칸을 갈 수 있는 경우도 있을 수 있으므로, max_step을 정할 수 있도록 했습니다.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct LatticeTopology{
    /// Maximum step size available in the system
    max_step : usize,
}

impl Default for LatticeTopology{
    fn default() -> Self{
        Self{
            max_step : 1usize,
        }
    }
}

impl LatticeTopology{
    /// Initializing [`LatticeTopology`](struct@LatticeTopology) with a maximum step size.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::topology::LatticeTopology;
    /// let sys = LatticeTopology::new(1);
    /// ```
    pub fn new(max_step : usize) -> Self{
        Self{
            max_step,
        }
    }
}

impl Topology<Cartessian1D<i32>> for LatticeTopology{
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
    /// # use moldybrody::system::topology::LatticeTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = LatticeTopology::new(1);
    /// let vt = Cartessian1D{coord : 0};
    ///
    /// let move1 = Cartessian1D{coord : 1};
    /// assert_eq!(sys.check_move(&vt, &move1), true);
    ///
    /// let move2 = Cartessian1D{coord : 2};
    /// assert_eq!(sys.check_move(&vt, &move2), false);
    /// ```
    ///
    /// # Panic
    ///
    /// point와 movement가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
    /// ```should_panic
    /// # use moldybrody::system::point::{Cartessian2D, CartessianND};
    /// # use moldybrody::system::topology::LatticeTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = LatticeTopology::new(1);
    /// let p = Cartessian2D{coord : [1; 2]};
    /// let move1 = CartessianND{coord : vec![0; 2]};
    /// // sys.check_move(&p, &move1);    // Cannot compile
    ///
    /// let p2 = CartessianND{coord : vec![0; 2]};
    /// let move2 = CartessianND{coord : vec![0, 1, 1]};
    /// sys.check_move(&p2, &move2);       // Panic!
    /// ```
    fn check_move(&self, _pos : &Cartessian1D<i32>, movement : &Cartessian1D<i32>) -> bool{
        let r = movement.taxi_norm();
        r <= self.max_step
    }
}


impl Topology<CartessianND<i32>> for LatticeTopology{
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
    /// # use moldybrody::system::topology::LatticeTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = LatticeTopology::new(1);
    /// let vt = CartessianND{coord : vec![0; 3]};
    ///
    /// let move1 = CartessianND{coord : vec![1, 0, 0]};
    /// assert_eq!(sys.check_move(&vt, &move1), true);
    ///
    /// let move2 = CartessianND{coord : vec![1, 1, 0]};
    /// assert_eq!(sys.check_move(&vt, &move2), false);
    /// ```
    ///
    /// # Panic
    ///
    /// point와 movement가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
    /// ```should_panic
    /// # use moldybrody::system::point::{Cartessian2D, CartessianND};
    /// # use moldybrody::system::topology::LatticeTopology;
    /// # use moldybrody::system::Topology;
    /// let sys = LatticeTopology::new(1);
    /// let p = Cartessian2D{coord : [1; 2]};
    /// let move1 = CartessianND{coord : vec![0; 2]};
    /// // sys.check_move(&p, &move1);    // Cannot compile
    ///
    /// let p2 = CartessianND{coord : vec![0; 2]};
    /// let move2 = CartessianND{coord : vec![0, 1, 1]};
    /// sys.check_move(&p2, &move2);       // Panic!
    /// ```
    fn check_move(&self, pos : &CartessianND<i32>, movement : &CartessianND<i32>) -> bool{
        if pos.dim() != movement.dim(){
            panic!("{}", ErrorCode::InvalidDimension);
        }
        let r = movement.taxi_norm();
        r <= self.max_step
    }
}

#[allow(unused_macros)]
macro_rules! impl_lattice_topology_nD{
    ($cartessian_name:ident, $dim:expr) =>{
        impl Topology<$cartessian_name<i32>> for LatticeTopology{

            doc_comment!{
                concat!(
                    "Check whether the movement is valid or not

point의 이동이 가능한 move인지 아닌지 여부를 확인해주는 함수이다.
시스템은 max_step이 정해져있는데,
이는 시뮬레이션 과정에서 입자가 움직일 수 있는 최대 변위를 제한한다.

# Examples

```
# use moldybrody::system::point::", stringify!($cartessian_name), ";
# use moldybrody::system::topology::LatticeTopology;
# use moldybrody::system::Topology;
let sys = LatticeTopology::new(1);
let vt = ", stringify!($cartessian_name), "{coord : [2; ", $dim, "]};
let mut move1 = ", stringify!($cartessian_name), "{coord : [0; ", $dim, "]};

move1.coord[0] = 1;   // length of movement is equal to 1
assert_eq!(sys.check_move(&vt, &move1), true);

move1.coord[1] = 1;       // length of movement is larger than 1
assert_eq!(sys.check_move(&vt, &move1), false);
```

# Panic

point와 movement가 서로 차원이 다르면 compile error가 생기거나 panic이 일어납니다.
```should_panic
# use moldybrody::system::point::{Cartessian2D, CartessianND};
# use moldybrody::system::topology::LatticeTopology;
# use moldybrody::system::Topology;
let sys = LatticeTopology::new(1);
let p = Cartessian2D{coord : [1; 2]};
let move1 = CartessianND{coord : vec![0; 2]};
// sys.check_move(&p, &move1);    // Cannot compile

let p2 = CartessianND{coord : vec![0; 2]};
let move2 = CartessianND{coord : vec![0; 3]};
sys.check_move(&p2, &move2);       // Panic!
```"
                ),
                fn check_move(&self, _pos : &$cartessian_name<i32>, movement : &$cartessian_name<i32>) -> bool{
                    let r = movement.taxi_norm();
                    r <= self.max_step
                }
            }
        }
    }
}

impl_lattice_topology_nD!(Cartessian2D, 2);
impl_lattice_topology_nD!(Cartessian3D, 3);
impl_lattice_topology_nD!(Cartessian4D, 4);

// ===========================================================================================
// ===========================================================================================
// ===========================================================================================

/// Network topology
///
/// Network의 topology는 node의 개수와 각 node들의 연결정보로 주어집니다.
/// 이때 연결정보는 연결 여부만 담는 경우가 있고,
/// edge에 weight가 부여된 경우가 있습니다.
/// 이를 포괄하기 위해 LatticeTopology의 adjoint matrix를 구성하는 type이 bool, i32, f64로 달라질 수 있도록 개발했습니다.
/// directed network의 경우에는 adjoint matrix가 symmetric하지 않을 수 있습니다.
/// 하지만 모든 경우에 adjoint matrix는 square matrix여야합니다. 그렇지 않으면 처음 정의에서부터 panic합니다. 
#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct NetworkTopology<T>{
    /// The number of nodes
    num_node : usize,
    /// Adjoint matrix,
    adjoint : Vec<Vec<T>>,
}

impl<T : Clone> NetworkTopology<T>{
    /// Initialization of Network toopology.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::topology::NetworkTopology;
    /// let adj_bool : Vec<Vec<bool>> = vec![vec![true, false], vec![true, true]];
    /// let net_bool = NetworkTopology::new(&adj_bool);
    ///
    /// let adj_int = vec![vec![1, 0], vec![1, 1]];
    /// let net_int = NetworkTopology::new(&adj_int);
    ///
    /// let adj_float = vec![vec![1.0, 0.0], vec![0.2, 0.5]];
    /// let net_float = NetworkTopology::new(&adj_float);
    /// ```
    ///
    /// # Panic
    /// Square matrix가 아닌 matrix가 들어오면 panic합니다. 
    /// ```should_panic
    /// # use moldybrody::system::topology::NetworkTopology;
    /// let adj_bool : Vec<Vec<bool>> = vec![vec![true, false], 
    ///                                      vec![true, true], 
    ///                                      vec![true, false]];
    /// let net_bool = NetworkTopology::new(&adj_bool);
    /// ```
    pub fn new(adjoint : &Vec<Vec<T>>) -> Self{
        // test square matrix
        let n = adjoint.len();
        for row in adjoint{
            if row.len() != n{
                panic!("{}", ErrorCode::InvalidArgumentInput);
            }
        }
        Self{
            num_node : adjoint.len(),
            adjoint : adjoint.clone(),
        }
    }
}

impl<T : Copy + Default + PartialEq> Topology<NodeIndex> for NetworkTopology<T>{
    /// Check whether the movement is valid or not
    ///
    /// point의 이동이 가능한 move인지 아닌지 여부를 확인해주는 함수이다.
    /// network에서는 목적지 node가 연결되어있는지 여부를 확인하는 방법으로 movement의 적합성을 확인합니다.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::NodeIndex;
    /// # use moldybrody::system::topology::NetworkTopology;
    /// # use moldybrody::system::Topology;
    /// let adj_bool : Vec<Vec<bool>> = vec![vec![true, false], vec![true, true]];
    /// let sys = NetworkTopology::new(&adj_bool);
    /// let v1 = NodeIndex{index : 0};
    /// let v2 = NodeIndex{index : 1};
    ///
    /// assert_eq!(sys.check_move(&v1, &v2), false);
    /// assert_eq!(sys.check_move(&v2, &v1), true);
    /// ```
    ///
    /// # Panic
    ///
    /// 출발점이나 도착점의 index가 시스템의 number of node를 넘어가면 panic합니다.
    /// ```should_panic
    /// # use moldybrody::system::point::NodeIndex;
    /// # use moldybrody::system::topology::NetworkTopology;
    /// # use moldybrody::system::Topology;
    /// let adj_bool : Vec<Vec<bool>> = vec![vec![true, false], vec![true, true]];
    /// let sys = NetworkTopology::new(&adj_bool);
    ///
    /// let v1 = NodeIndex{index : 0};
    /// let v2 = NodeIndex{index : 3};
    /// sys.check_move(&v1, &v2);           // Panic!
    /// sys.check_move(&v2, &v1);           // Panic!
    /// ```
    fn check_move(&self, depart : &NodeIndex, arrival : &NodeIndex) -> bool{
        let connect = self.adjoint[depart.index][arrival.index];
        !(connect == T::default())
    }
}
