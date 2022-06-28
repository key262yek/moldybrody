//! Rust package for Molecular dynamics and Brownian Dynamics simulation
//!
//! 이 package는 Molecular Dynamics(MD)와 Brownian Dynamics(BD) simulation code에 필수적인 요소들의 추상적 구조를 제공합니다.
//! - [System](crate::system)는 입자들이 돌아다니는 공간의 다양한 topology 구조와 그에 필요한 함수들을 제공합니다.
//! - [Boundary Condition](crate::boundary)은 주어진 공간의 다양한 경계조건을 구현할 수 있도록 합니다.
//! - [Force computation](crate::boundary)은 입자들이 받게 되는 상호작용을 계산하는 방식을 추상화한 module로, Global potential, Two body interaction 뿐 아니라 다양한 random force까지 포괄할 수 있습니다.
//! - [Approximation method](crate::approx)는 앞서 계산한 Force를 이용해 입자들의 위치 및 운동 정보를 갱신하는 방법론들을 제공합니다.
//! - [Data Analysis](crate::analysis)는 시뮬레이션을 통해 얻은 결과를 분석하는 도구들을 제공합니다.
//!
//! 아래의 example을 통해 시뮬레이션 구성을 어떻게 간소화했는지 확인해보세요.
//! ```rust
//!   let x = 10;
//! ```

extern crate doc_comment;
extern crate ndarray;
extern crate num_traits;

pub mod approx;
pub mod boundary;
pub mod force;
pub mod state;
pub mod vector;
// pub mod analysis;
pub mod iterator;
pub mod rng;

pub(crate) mod error;
pub mod prelude;

//
// Todo
//
// Common
// Argument information auto-complete
// Argument parsing
//
// System
// Continuous : Rectangular, Circular
// Point & Boundary로 정의됨.
// Boundary condition check
// Periodic
// Reflective
// Robin
// etc
// Point
// Associated with System
// Vector operation : Add, Scalar Multiplication, Sub, Dot, cross
// Boundary
// 경계면은 대부분 함수로 정의됨.
// Agent
// Single-ptl : Global Interaction, Brownian, Active motion, Levy walk
// Double-ptl : Central forces
// Move
// Approximation
// First order approximation : Newtonian / Stochastic
// Second .. : Newtonian / Stochastic
// Forces to Displacement
// Analysis
