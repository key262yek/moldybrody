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
pub trait Point{

}

/// Neighbor of a point
///
/// Topology를 구성하는 요소 중 이웃에 해당하는 trait.
/// 특정 점이 이웃에 속하는지 여부를 확인할 수 있어야 합니다.
pub trait Neighbor<'a>{

}


/// Topology on the system
///
/// Point와 Neighbor trait을 구성요소로 하여 시스템의 Topology를 정의해주는 trait입니다.
/// 점과 이웃의 관계는 비로소 Topology 안에서 정의됩니다.
pub trait Topology<'a, P, N>
    where P : Point,
          N : Neighbor<'a>{
    /// Return a neighbor of the point
    ///
    /// 점의 이웃을 알려주는 함수입니다.
    /// 시뮬레이션에서 이웃의 보다 실천적 의미는 한 번에 갈 수 있는 점들의 집합을 의미합니다.
    /// 
    /// # Examples 
    /// 
    /// ```
    /// use moldybrody::system::point::Cartessian3D;
    /// use moldybrody::system::neighbor::OpenBall3D;
    /// use moldybrody::system::topology::ContinuousTopology;
    /// use moldybrody::system::Topology;
    ///
    /// let sys = ContinuousTopology::new(0.1);
    /// let v = Cartessian3D{coord : [0.0; 3]};
    /// assert_eq!(sys.neighbor(&v), OpenBall3D::new(&v));
    /// ```
    fn neighbor<'b : 'a>(&self, pos : &'b P) -> N;

    /// Check whether a point is in the neighbor or not
    ///
    /// point가 neighbor에 속하는지 여부를 확인하는 함수입니다.
    /// particle이 특정 위치로 이동할 수 있는지 여부를 test할 때 쓰입니다.
    /// 
    /// # Examples 
    /// 
    /// ```
    /// use moldybrody::system::point::Cartessian2D;
    /// use moldybrody::system::neighbor::OpenBall2D;
    /// use moldybrody::system::topology::ContinuousTopology;
    /// use moldybrody::system::Topology;
    ///
    /// let sys = ContinuousTopology::new(0.1);
    /// let vt = Cartessian2D{coord : [0.0, 0.05]};
    ///
    /// let v1 = Cartessian2D{coord : [0.0, 0.0]};
    /// let n1 = OpenBall2D::new(&v1);
    /// assert_eq!(sys.inclusion(&vt, &n1), true);
    ///
    /// let v2 = Cartessian2D{coord : [0.0, 1.0]};
    /// let n2 = OpenBall2D::new(&v2);
    /// assert_eq!(sys.inclusion(&vt, &n2), false);
    /// ```
    fn inclusion(&self, pos : &P, neigh : &N) -> bool;
}


pub mod point;
pub mod neighbor;
pub mod topology;

