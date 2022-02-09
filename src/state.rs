//! Implement general form of state of particles
//!
//! 시뮬레이션에서 계속 갱신하는 입자의 상태는 시스템에 따라 그 정보의 양이 결정됩니다.
//! 일반적인 시스템에서는 입자의 위치와 속도를 모두 알아야하지만,
//! overdamped regime에서는 입자의 위치를 아는 것만으로도 속도를 유추할 수 있습니다.
//! 따라서 이 module에서는 일반적인 상황에서의 상태와 overdamped regime에서의 상태를 각기 정의해놓고자 합니다.

use crate::vector::Vector;

pub trait State{
    type Movement : Vector;
    type Position : Vector;

    fn pos(&self) -> &Self::Position;
    fn pos_mut(&mut self) -> &mut Self::Position;

    fn disp<'a>(&self, movement : &'a Self::Movement) -> &'a Self::Position;

    fn renew_state(&mut self, movement : &Self::Movement);
}


pub trait HasVelocity : State{

    fn vel(&self) -> &<Self as State>::Position;
    fn vel_mut(&mut self) -> &mut <Self as State>::Position;
}

pub trait Mass<T>{
    fn mass(&self) -> T;
}

pub trait Charge<T>{
    fn charge(&self) -> T;
}

// pub trait GeneralState : State<T, N>{
//     type Position : Vector;
//     type Velocity : Vector;

//     fn pos(&self) -> &Self::Position;
//     fn pos_mut(&mut self) -> &mut Self::Position;

//     fn vel(&self) -> &Self::Velocity;
//     fn vel_mut(&mut self) -> &mut Self::Velocity;
// }


// pub trait OverdampedState : State{
//     type Position : Vector;

//     fn pos(&self) -> &Self::Position;
//     fn pos_mut(&mut self) -> &mut Self::Position;
// }

#[cfg(test)]
mod test {
    use crate::vector::product::Norm;
    use super::*;
    use crate::prelude::{Mass, Charge, State};
    use crate::vector::{Cartessian2D};

    #[test]
    fn test_derive_macro(){

        #[derive(Mass, Charge)]
        struct TestState{
            mass : f64,
            charge : f64,
            #[allow(dead_code)]
            pos : Cartessian2D<f64>,
        }

        let test = TestState{mass : 0f64, charge : 0f64, pos : Cartessian2D::new([0f64, 0f64])};
        assert!(test.mass() < 1e-30);
        assert!(test.charge() < 1e-30);


        #[derive(State)]
        struct TestState2{
            mass : f64,
            charge : f64,
            pos : Cartessian2D<f64>,
        }
        let mut test = TestState2{mass : 0f64, charge : 0f64, pos : Cartessian2D::new([0f64, 0f64])};
        assert!(test.mass() < 1e-30);
        assert!(test.charge() < 1e-30);
        assert!(test.pos().norm_l2() < 1e-30);
        let movement = Cartessian2D::new([1f64, 2f64]);
        test.renew_state(&movement);
        assert!((test.pos().norm_l2() - 5f64.sqrt()).abs() < 1e-30);

        #[derive(State)]
        struct TestState3{
            mass : f64,
            charge : f64,
            pos : Cartessian2D<f64>,
            vel : Cartessian2D<f64>,
        }
        let mut test = TestState3{mass : 0f64,
                    charge : 0f64,
                    pos : Cartessian2D::new([0f64, 0f64]),
                    vel : Cartessian2D::new([0f64, 0f64])};
        assert!(test.mass() < 1e-30);
        assert!(test.charge() < 1e-30);
        assert!(test.pos().norm_l2() < 1e-30);
        assert!(test.vel().norm_l2() < 1e-30);
        let movement = (Cartessian2D::new([1f64, 2f64]), Cartessian2D::new([0f64, 1f64]));
        test.renew_state(&movement);
        assert!((test.pos().norm_l2() - 5f64.sqrt()).abs() < 1e-30);
        assert!((test.vel().norm_l2() - 1f64).abs() < 1e-30);
    }


}

