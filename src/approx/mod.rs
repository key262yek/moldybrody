//! Approximation method to expect time variation of state
//!
//! Molecular dynamics도, Stochastic dynamics도 모두 여러 approximation method가 존재합니다.
//! 각 dynamics에서, general case냐 overdamped case냐에 따라, 그리고 order에 따라 달라지는 approximation method를 정리하였습니다.

use crate::prelude::*;

/// Indicate a struct as state of particle
///
/// 어떤 strcut는 각 함수들에서 particle의 state 역할을 할 것입니다.
/// state struct란 것을 지칭하기 위한 trait입니다.
pub trait State {}

/// Approximation method
///
/// System에 따라 정해지는 State와 Force, Time Iterator를 인자로 해서 State를 renewal하는 approximation method를 지칭합니다.
/// 같은 시스템에서도 approximation order에 따라 그 형태가 달라질 수 있습니다.
pub trait Approximation<S, P, T>
    where Self : Sized,
          S : State,
          P : Point{
    /// Renew state by using approximation
    fn renew(&self, state : &mut S, force : &P, dt : T);
}

/// Trait for iterate time of system
///
/// Simulation에서는 time을 계속 iterate할 필요성이 있습니다.
/// 어떤 경우는 time step이 계속 일정할 수도 있고,
/// 초반의 time step은 지수함수적으로 증가하는 경우도 있을 수 있습니다.
/// 또 입자들이 Run and Tumble하는 경우에는, 다음 tumble까지 직진운동할 것이기 때문에 직진운동 사이에서는 시간을 크게 뛸 수 있습니다.
/// 이런 여러 성격의 iterator를 모두 어우르기 위해 trait으로 선언했습니다.
pub trait TimeIterator<T>
    where Self : Sized{
    /// Return a current time restore in timeiterator
    ///
    /// Iterator에 저장되어 있는 현재 시간을 출력하는 함수입니다.
    ///
    /// # Examples
    /// ```
    /// # use moldybrody::approx::time::ConstStep;
    /// # use moldybrody::approx::TimeIterator;
    /// let const_step = ConstStep::new(1e-3).unwrap();
    /// assert_eq!(const_step.current_time(), const_step.current);
    /// ```
    fn current_time(&self) -> T;

    /// Return a time step
    ///
    /// 현재 시점에 Iterator의 time step을 출력하는 함수
    ///
    /// # Examples
    /// ```
    /// # use moldybrody::approx::time::ConstStep;
    /// # use moldybrody::approx::TimeIterator;
    /// let const_step = ConstStep::new(1e-3).unwrap();
    /// assert_eq!(const_step.dt, 1e-3);
    /// ```
    fn dt(&self) -> T;

    /// Change upper limit of time
    ///
    /// Iterator에는 time의 upper limit을 정의할 필요가 있는데,
    /// 그 maximum time을 바꾸는 함수입니다.
    ///
    /// # Examples
    /// ```
    /// # use moldybrody::approx::time::ConstStep;
    /// # use moldybrody::approx::TimeIterator;
    /// let mut const_step = ConstStep::new(1e-3).unwrap();
    ///
    /// const_step.set_tmax(1f64);
    /// assert_eq!(const_step.tmax, 1f64);
    ///
    /// const_step.set_tmax(10f64);
    /// assert_eq!(const_step.tmax, 10f64);
    /// ```
    fn set_tmax(&mut self, tmax : T) -> Result<()>;

    /// Return iterator of both a current time and time step simultaneously
    ///
    /// TimeIterator는 current time의 iterator입니다.
    /// 하지만 simulation에는 current time 뿐 아니라 time step 역시 필요합니다.
    /// 이 둘을 동시에 iterate해주는 iterator가 TimeDiffIterator입니다.
    ///
    /// # Examples
    /// ```ignore
    /// # use moldybrody::approx::time::ConstStep;
    /// # use moldybrody::approx::TimeIterator;
    /// let mut const_step = ConstStep::new(1e-3).unwrap();
    /// for time in const_step{
    ///    // code with only time
    /// }
    ///
    /// const_step.reset();
    /// for (time, dt) in const_step.into_diff(){
    ///    // code with time and dt
    /// }
    /// ```

    fn into_diff(&self) -> TimeDiffIterator<Self>;
}

/// Iterator of both current time and time step simultaneously
///
/// 현재 시점과 시간간격을 동시에 iterate해주는 iterator입니다.
///
/// # Examples
/// ```
/// # use moldybrody::approx::time::ConstStep;
/// # use moldybrody::approx::TimeIterator;
/// let const_step = ConstStep::new(1e-3).unwrap();
/// let mut const_diff = const_step.into_diff();
/// assert_eq!(const_diff.next(), Some((0f64, 1e-3)));
/// assert_eq!(const_diff.next(), Some((1e-3, 1e-3)));
/// assert_eq!(const_diff.next(), Some((2e-3, 1e-3)));
/// ```
///
/// ```ignore
/// # use moldybrody::approx::time::ConstStep;
/// # use moldybrody::approx::TimeIterator;
/// let mut const_step = ConstStep::new(1e-3).unwrap();
/// for (time, dt) in const_step.into_diff(){
///     // something with time and dt
/// }
/// ```
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct TimeDiffIterator<T>{
    pub timeiter : T,
}

impl<T : Iterator + TimeIterator<f64>> Iterator for TimeDiffIterator<T>{
    type Item = (T::Item, f64);

    fn next(&mut self) -> Option<Self::Item>{
        match self.timeiter.next(){
            Some(time) => Some((time, self.timeiter.dt())),
            None => None,
        }
    }
}

pub mod state;
pub mod time;
pub mod molecular_dynamics;
pub mod brownian_dynamics;
