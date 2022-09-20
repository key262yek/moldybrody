//! Implement general form of state of particles
//!
//! 시뮬레이션에서 계속 갱신하는 입자의 상태는 시스템에 따라 그 정보의 양이 결정됩니다.
//! 일반적인 시스템에서는 입자의 위치와 속도를 모두 알아야하지만,
//! overdamped regime에서는 입자의 위치를 아는 것만으로도 속도를 유추할 수 있습니다.
//! 따라서 이 module에서는 일반적인 상황에서의 상태와 overdamped regime에서의 상태를 각기 정의해놓고자 합니다.

use crate::vector::Vector;

pub trait State {
    type Movement: Vector;
    type Position: Vector;

    fn pos(&self) -> &Self::Position;
    fn pos_mut(&mut self) -> &mut Self::Position;

    fn disp<'a>(&self, movement: &'a Self::Movement) -> &'a Self::Position;

    fn renew_state(&mut self, movement: &Self::Movement);

    fn renew_with_constant(&mut self, movement : &Self::Movement, c : <Self::Position as Vector>::Item);
}

pub trait HasVelocity: State {
    fn vel(&self) -> &<Self as State>::Position;
    fn vel_mut(&mut self) -> &mut <Self as State>::Position;
}

pub trait Mass<T>: State {
    fn mass(&self) -> T;
}

pub trait Charge<T>: State {
    fn charge(&self) -> T;
}

pub trait DiffusionFloat<T>: State {
    fn diff_const(&self) -> T;
}

pub trait DiffusionInt<T>: State {
    fn diff_time(&self) -> T;
}

pub trait Orientation<P>: State
where
    P: Vector,
{
    fn orientation(&self) -> P;
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
    use super::*;
    use crate::prelude::State;
    use crate::vector::product::Norm;
    use crate::vector::basic::Map;
    use crate::vector::Cartessian2D;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_derive_macro() {
        #[derive(State)]
        struct TestState {
            mass: f64,
            charge: f64,
            #[allow(dead_code)]
            pos: Cartessian2D<f64>,
        }

        let test = TestState {
            mass: 0f64,
            charge: 0f64,
            pos: Cartessian2D::new([0f64, 0f64]),
        };
        assert_abs_diff_eq!(test.mass(), 0f64, epsilon = 1e-30);
        assert_abs_diff_eq!(test.charge(), 0f64, epsilon = 1e-30);

        #[derive(State)]
        struct TestState2 {
            mass: f64,
            charge: f64,
            pos: Cartessian2D<f64>,
        }
        let mut test = TestState2 {
            mass: 0f64,
            charge: 0f64,
            pos: Cartessian2D::new([0f64, 0f64]),
        };
        assert_abs_diff_eq!(test.mass(), 0f64, epsilon = 1e-30);
        assert_abs_diff_eq!(test.charge(), 0f64, epsilon = 1e-30);
        assert_abs_diff_eq!(test.pos().norm_l2(), 0f64, epsilon = 1e-30);
        let movement = Cartessian2D::new([1f64, 2f64]);
        test.renew_state(&movement);
        assert_abs_diff_eq!(test.pos(), &movement, epsilon = 1e-30);

        #[derive(State)]
        struct TestState3 {
            mass: f64,
            charge: f64,
            pos: Cartessian2D<f64>,
            vel: Cartessian2D<f64>,
        }
        let mut test = TestState3 {
            mass: 0f64,
            charge: 0f64,
            pos: Cartessian2D::new([0f64, 0f64]),
            vel: Cartessian2D::new([0f64, 0f64]),
        };
        assert_abs_diff_eq!(test.mass(), 0f64, epsilon = 1e-30);
        assert_abs_diff_eq!(test.charge(), 0f64, epsilon = 1e-30);
        assert_abs_diff_eq!(test.pos().norm_l2(), 0f64, epsilon = 1e-30);
        assert_abs_diff_eq!(test.vel().norm_l2(), 0f64, epsilon = 1e-30);
        let movement = (
            Cartessian2D::new([1f64, 2f64]),
            Cartessian2D::new([0f64, 1f64]),
        );
        test.renew_state(&movement);
        assert_abs_diff_eq!(test.pos(), &movement.0, epsilon = 1e-30);
        assert_abs_diff_eq!(test.vel(), &movement.1, epsilon = 1e-30);

        test.renew_with_constant(&movement, -2f64);
        assert_abs_diff_eq!(test.pos(), &(-movement.0), epsilon = 1e-30);
        assert_abs_diff_eq!(test.vel(), &(-movement.1), epsilon = 1e-30);

        // #[derive(State)]
        // struct TestState4 {
        //     diff_const: f64,
        //     pos: Cartessian2D<f64>,
        // }
        // let test = TestState4 {
        //     diff_const: 1f64,
        //     pos: Cartessian2D::new([0f64; 2]),
        // };
        // assert_abs_diff_eq!(test.diff_const(), 1f64, epsilon = 1e-3);

        // #[derive(State)]
        // struct TestState5 {
        //     diff_const: [f64; 3],
        //     pos: Cartessian2D<f64>,
        // }

        // let test = TestState5 {
        //     diff_const: [1f64, 2f64, 3f64],
        //     pos: Cartessian2D::new([0f64; 2]),
        // };
        // assert_eq!(test.diff_const(), [1f64, 2f64, 3f64]);
    }
}
