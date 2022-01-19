//! Provide general structure of force and potential
//!
//! 입자들은 시스템에서 아래와 같은 여러 상호작용을 받을 수 있습니다.
//! - Global Potential : 시스템 전체에 존재하는 force field
//! - Bimolecular Interaction : 두 입자가 주고받는 상호작용
//! - Random Force : 주변 환경에 의해 받게 되는 random force
//!
//! Force, Potential trait을 정의해두어 훗날 변수들을 구분하기 위해 사용할 것이며,
//! GlobalPotential, BimolecularInteraction, RandomForce trait을 정의하여 각각의 경우에 필요로 하는 함수들을 정의하였습니다.

use crate::prelude::*;

/// Global potential, Single particle interaction
pub trait GlobalPotential<P>
    where P : Point{

    /// return value of potential at certain state
    fn potential(&self, state : &P) -> f64;


    /// change value of force to force vector of given state
    fn force_to(&self, state : &P, force : &mut P);
}

/// Trait indicates bimolecular interaction
pub trait BimolecularInteraction<P>
    where P : Point{

    /// return value of potential at certain state
    fn potential(&self, disp : &P) -> f64;

    /// change value of force to force vector of given state
    fn force_to(&self, disp : &P, force : &mut P);
}

/// Trait indicates random force
pub trait RandomForce<P : Point>{
    /// change value of force to random force vector of given state
    fn force_to(&self, rng : &mut Pcg64, force : &mut P);
}

pub mod global;
pub mod bimolecular;
pub mod random;
