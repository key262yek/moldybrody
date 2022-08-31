//! Approximation method to expect time variation of state
//!
//! Molecular dynamics도, Stochastic dynamics도 모두 여러 approximation method가 존재합니다.
//! 각 dynamics에서, general case냐 overdamped case냐에 따라, 그리고 order에 따라 달라지는 approximation method를 정리하였습니다.

// use crate::prelude::*;

use crate::error::Error;
// use crate::force::Bimolecular;
// use crate::force::Global;
use crate::state::State;
use crate::vector::Scalar;
use crate::vector::Vector;
// use std::ops::{Index, IndexMut};
use approx::AbsDiffEq;

pub struct ButcherTableau<T>{
    aij : Vec<Vec<T>>,
    bi : Vec<T>,
    ci : Vec<T>,
}

enum NewtonButcherTableauBuilder<T>{
    Euler,
    MidPoint,
    Ssprk3,
    Classic,
    Heun(usize),
    Ralston(usize),
    Generic(usize, T),
}

impl<T> NewtonButcherTableauBuilder<T>{
    fn new(name : Option<&str>, order : Option<usize>, alpha : Option<T>) -> Self{
        match name{
            Some("Euler") => NewtonButcherTableauBuilder::Euler,
            Some("MidPoint") => NewtonButcherTableauBuilder::MidPoint,
            Some("SSPRK3") => NewtonButcherTableauBuilder::Ssprk3,
            Some("Classic") => NewtonButcherTableauBuilder::Classic,
            Some("Heun") => {
                match order {
                    Some(v) => NewtonButcherTableauBuilder::Heun(v),
                    None => panic!("Heun method requires order attribute."),
                }
            },
            Some("Ralston") => {
                match order {
                    Some(v) => NewtonButcherTableauBuilder::Ralston(v),
                    None => panic!("Ralston method requires order attribute."),
                }
            },
            Some("Generic") => {
                match (order, alpha) {
                    (Some(v), Some(a)) => NewtonButcherTableauBuilder::Generic(v, a),
                    (Some(_), None) => panic!("Generic Runge-Kutta method requires parameter alpha in [0, 1]"),
                    (None, _) => panic!("Generic Runge-Kutta method requires order"),
                }
            }
            None | Some(_) => panic!("Butcher Tableau builder requires valid name of tableau. ex) Euler, Midpoint, SSPRK3, Classic, Heun, Ralston, Generic"),
        }
    }
}

macro_rules! impl_bt_builder{
    ($ty : ident) => {
        impl NewtonButcherTableauBuilder<$ty> {
            fn build(self) -> ButcherTableau<$ty>{
                match self{
                    NewtonButcherTableauBuilder::Euler => {
                        ButcherTableau{
                            aij : vec![vec![0.]],
                            bi : vec![1.],
                            ci : vec![0.],
                        }
                    },
                    NewtonButcherTableauBuilder::MidPoint => {
                        ButcherTableau{
                            aij : vec![vec![0., 0.], vec![0.5, 0.]],
                            bi : vec![0., 1.],
                            ci : vec![0., 0.5],
                        }
                    },
                    NewtonButcherTableauBuilder::Ssprk3 => {
                        let t = 1. / 6.;
                        ButcherTableau{
                            aij : vec![vec![0., 0., 0.], vec![1., 0., 0.], vec![0.25, 0.25, 0.]],
                            bi : vec![t, t, 4. * t],
                            ci : vec![0., 1., 0.5],
                        }    
                    },
                    NewtonButcherTableauBuilder::Classic => {
                        let t = 1. / 6.;
                        ButcherTableau{
                            aij : vec![vec![0., 0., 0., 0.], vec![0.5, 0., 0., 0.], vec![0., 0.5, 0., 0.], vec![0., 0., 1., 0.]],
                            bi : vec![t, 2. * t, 2. * t, t],
                            ci : vec![0., 0.5, 0.5, 1.],
                        } 
                    },
                    NewtonButcherTableauBuilder::Heun(order) => {
                        match order {
                            2 => {
                                ButcherTableau{
                                    aij : vec![vec![0., 0.], vec![1., 0.]],
                                    bi : vec![0.5, 0.5],
                                    ci : vec![0., 1.],
                                }
                            },
                            3 => {
                                let t = 1. / 3.;
                                ButcherTableau{
                                    aij : vec![vec![0., 0., 0.], vec![t, 0., 0.], vec![0., 2. * t, 0.]],
                                    bi : vec![0.25, 0., 0.75],
                                    ci : vec![0., t, 2. * t],
                                }
                            }, 
                            _ => {
                                panic!("Heun methods exist only for 2nd and 3rd order");
                            }
                        }
                    },
                    NewtonButcherTableauBuilder::Ralston(order) => {
                        match order {
                            2 => {
                                ButcherTableau{
                                    aij : vec![vec![0., 0.], vec![2. / 3., 0.]],
                                    bi : vec![0.25, 0.75],
                                    ci : vec![0., 2. / 3.],
                                }
                            },
                            3 => {
                                let t = 1. / 9.;
                                ButcherTableau{
                                    aij : vec![vec![0., 0., 0.], vec![0.5, 0., 0.], vec![0., 0.75, 0.]],
                                    bi : vec![2. * t, 3. * t, 4. * t],
                                    ci : vec![0., 0.5, 0.75],
                                }    
                            }, 
                            _ => {
                                panic!("Raltson methods exist only for 2nd and 3rd order");
                            }
                        }
                    },
                    NewtonButcherTableauBuilder::Generic(order, alpha) => {
                        match order {
                            2 => {
                                let t = 1. / (2. * alpha);
                                ButcherTableau{
                                    aij : vec![vec![0., 0.], vec![alpha, 0.]],
                                    bi : vec![1. - t, t],
                                    ci : vec![0., alpha],
                                }
                            },
                            3 => {
                                if alpha.abs_diff_eq(&0., 1e-5) || alpha.abs_diff_eq(&(2. / 3.), 1e-5) || alpha.abs_diff_eq(&1., 1e-5){
                                    panic!("Generic 3rd order Runge-Kutta method cannot be made with alpha = 0, 2/3, 1");
                                }
                                let t = (1. - alpha) / (3. * alpha - 2.);
                                ButcherTableau{
                                    aij : vec![vec![0., 0., 0.], vec![alpha, 0., 0.], vec![1. + t / alpha, - t / alpha, 0.]],
                                    bi : vec![0.5 - 1. / (6. * alpha), 1. / (6. * alpha * (1. - alpha)), - t / 6.],
                                    ci : vec![0., alpha, 1.],
                                }
                            }, 
                            _ => {
                                panic!("Generic Runge-Kutta tableau exist only for 2nd and 3rd order");
                            }
                        }
                    }
                }
            }
        }
    }
}
impl_bt_builder!(f32);
impl_bt_builder!(f64);



/// Approximation method
///
/// System에 따라 정해지는 State와 Force, Time Iterator를 인자로 해서 State를 renewal하는 approximation method를 지칭합니다.
/// 같은 시스템에서도 approximation order에 따라 그 형태가 달라질 수 있습니다.
pub struct NewtonRKIntegratorBuilder<T>{
    name_tableau : Option<&'static str>,
    order_tableau : Option<usize>,
    alpha_tableau : Option<T>,
}
impl<T : Copy> NewtonRKIntegratorBuilder<T>{
    pub fn new() -> Self{
        Self{
            name_tableau : None,
            order_tableau : None,
            alpha_tableau : None,
        }
    }

    pub fn name_tableau(mut self, name : &'static str) -> Self{
        self.name_tableau = Some(name);
        self
    }

    pub fn order_tableau(mut self, order : usize) -> Self{
        self.order_tableau = Some(order);
        self
    }

    pub fn alpha_tableau(mut self, alpha : T) -> Self{
        self.alpha_tableau = Some(alpha);
        self
    }

    pub fn build() -> NewtonRKIntegrator{
        unimplemented!();
    }
}

pub struct NewtonRKIntegrator(Box<dyn NewtonRKIntegratorTrait>);
trait NewtonRKIntegratorTrait{
    fn iterate(&self){
        unimplemented!();
    }
}


pub trait ApproxOverdampedLangevin<'a, P, T>: State<Movement = P, Position = P>
where
    P: Vector<Item = T>,
    T: Scalar,
{
    /// Pure diffusion
    fn euler_pure_diffusion(&'a self, random_force: &'a P, dt: T) -> P;

    fn euler_with_drift(&'a self, force: &'a P, random_force: &'a P, dt: T) -> P;
}

/// Trait for iterate time of system
///
/// Simulation에서는 time을 계속 iterate할 필요성이 있습니다.
/// 어떤 경우는 time step이 계속 일정할 수도 있고,
/// 초반의 time step은 지수함수적으로 증가하는 경우도 있을 수 있습니다.
/// 또 입자들이 Run and Tumble하는 경우에는, 다음 tumble까지 직진운동할 것이기 때문에 직진운동 사이에서는 시간을 크게 뛸 수 있습니다.
/// 이런 여러 성격의 iterator를 모두 어우르기 위해 trait으로 선언했습니다.
pub trait TimeIterator<T>
where
    Self: Sized,
{
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
    fn set_tmax(&mut self, tmax: T) -> Result<(), Error>;

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
pub struct TimeDiffIterator<T> {
    pub timeiter: T,
}

impl<T: Iterator + TimeIterator<f64>> Iterator for TimeDiffIterator<T> {
    type Item = (T::Item, f64);

    fn next(&mut self) -> Option<Self::Item> {
        match self.timeiter.next() {
            Some(time) => Some((time, self.timeiter.dt())),
            None => None,
        }
    }
}

pub mod langevin;
pub mod newton;
pub mod time;
// pub mod brownian_dynamics;
