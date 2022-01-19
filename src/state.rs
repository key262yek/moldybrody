//! Implement general form of state of particles
//!
//! 시뮬레이션에서 계속 갱신하는 입자의 상태는 시스템에 따라 그 정보의 양이 결정됩니다.
//! 일반적인 시스템에서는 입자의 위치와 속도를 모두 알아야하지만,
//! overdamped regime에서는 입자의 위치를 아는 것만으로도 속도를 유추할 수 있습니다.
//! 따라서 이 module에서는 일반적인 상황에서의 상태와 overdamped regime에서의 상태를 각기 정의해놓고자 합니다.

// use crate::vector::Vector;

// pub trait State<T, const N : usize> {}

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


