//! Construct a basic structure of system topology, and provide widely used topologies.
//!
//! 입자가 움직이는 공간은 공간 속 '점'과 그 점에서 옮겨갈 수 있는 '이웃'으로 정의됩니다.
//! 즉, 공간의 Topology가 곧 공간의 특성을 좌우합니다.
//! system module에서는 이 지점에 착안해 공간의 topology를 다양하게 정의했습니다.
//! 크게 세 개의 topology가 있습니다.
//! - Continuous system : open ball이 neighbor로 주어짐
//! - Lattice : connected node가 neighbor
//! - Complex network : connected node가 neighbor

/// Point in the system
///
/// Topology를 구성하는 요소 중 점에 해당하는 trait.
/// 각 점에 대응되는 이웃을 알 수 있어야합니다.
pub trait Point{}


/// Topology on the system
///
/// Point와 Neighbor trait을 구성요소로 하여 시스템의 Topology를 정의해주는 trait입니다.
/// 점과 이웃의 관계는 비로소 Topology 안에서 정의됩니다.
pub trait Topology<P>
    where P : Point{
    /// Check whether a movement is valid or not.
    ///
    /// movement나 목적지가 주어졌을 때, 해당 movement가 가능한 move인지 아닌지에 대해 확인하는 함수입니다.
    /// 
    /// # Examples 
    /// 
    /// ```
    /// use moldybrody::system::point::Cartessian2D;
    /// use moldybrody::system::topology::ContinuousTopology;
    /// use moldybrody::system::Topology;
    ///
    /// let sys = ContinuousTopology::new(0.1);
    /// let vt = Cartessian2D{coord : [0.0; 2]};
    ///
    /// let move1 = Cartessian2D{coord : [0.05, 0.0]};
    /// assert_eq!(sys.inclusion(&vt, &move1), true);
    ///
    /// let move2 = Cartessian2D{coord : [0.0, 1.1]};
    /// assert_eq!(sys.inclusion(&vt, &move2), false);
    /// ```
    fn check_move(&self, pos : &P, movement : &P) -> bool;
}

pub mod point;
pub mod topology;

