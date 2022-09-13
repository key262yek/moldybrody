//! Provide general structure of force and potential
//!
//! 입자들은 시스템에서 아래와 같은 여러 상호작용을 받을 수 있습니다.
//! - Global Potential : 시스템 전체에 존재하는 force field
//! - Bimolecular Interaction : 두 입자가 주고받는 상호작용
//! - Random Force : 주변 환경에 의해 받게 되는 random force
//!
//! Force, Potential trait을 정의해두어 훗날 변수들을 구분하기 위해 사용할 것이며,
//! GlobalPotential, BimolecularInteraction, RandomForce trait을 정의하여 각각의 경우에 필요로 하는 함수들을 정의하였습니다.

use crate::state::State;
use crate::vector::Vector;
use rand_pcg::Pcg64;
use std::fmt::Debug;

/// Global potential, Single particle interaction
pub trait Global<'a, S: State> : Debug{
    type Force: Vector;
    type Potential;

    /// return value of potential at certain state
    fn potential(&'a self, state: &'a S) -> Self::Potential;

    fn force(&'a self, state: &'a S) -> Self::Force;

    /// change value of force to force vector of given state
    fn force_to(&'a self, state: &'a S, force: &'a mut Self::Force);

    /// Add force to given vector
    fn force_add_to(&'a self, state: &'a S, force: &'a mut Self::Force);
}

/// Trait indicates bimolecular interaction
pub trait Bimolecular<'a, S: State> : Debug {
    type Force: Vector;
    type Potential;

    /// return value of potential at certain state
    fn potential(&'a self, state: &'a S, other: &'a S) -> Self::Potential;

    fn force(&'a self, state: &'a S, other: &'a S) -> Self::Force;

    /// change value of force to force vector of given state
    fn force_to(&'a self, state: &'a S, other: &'a S, force: &'a mut Self::Force);

    /// Add force to given vector
    fn force_add_to(&'a self, state: &'a S, other: &'a S, force: &'a mut Self::Force);
}

/// Trait indicates random force
pub trait RandomForce<'a, S: State> : Debug{
    type Force: Vector;

    fn force(&self, rng: &mut Pcg64) -> Self::Force;

    /// change value of force to random force vector of given state
    fn force_to(&self, rng: &mut Pcg64, force: &mut Self::Force);

    /// Add force to given vector
    fn force_add_to(&self, rng: &mut Pcg64, force: &mut Self::Force);
}

pub mod bimolecular;
pub mod global;
pub mod random;
