//! Implement widely used time iterators
//!
//! 많은 영역에서 사용되고 있는 time iterator들을 구현해놓았습니다.
//! - Constant Step : 일정한 time step으로 시뮬레이션할 때 사용할 수 있습니다.
//! - Exponential Step : 시스템 초반에는 짧은 time step으로 시뮬레이션하고, 그 이후부터는 constant step으로 시뮬레이션 할 때 사용할 수 있습니다.
//! - Nearest Time : Run and Tumble 방식으로 움직이는 입자들의 경우, Run 과정에서는 굳이 짧은 시간 단위로 끊지 않아도 오차 없이 입자들의 운동을 기술할 수 있습니다. 따라서 우리는 각 입자의 Tumble time을 ordered list에 보관하고 가장 근접한 시점으로 바로 이동하는 방식을 채택할 수 있습니다.

// use crate::prelude::*;
use crate::approx::TimeDiffIterator;
use crate::approx::TimeIterator;
use crate::error::{Error, ErrorCode};
use serde::{Deserialize, Serialize};

/// Constant step time iterator
///
/// 시스템 시간을 일정한 time step으로 변화시킬 때 사용하는 time iterator입니다.
///
/// # Examples
/// ```
/// # use moldybrody::approx::time::ConstStep;
/// # use moldybrody::approx::TimeIterator;
/// let mut const_step = ConstStep::new(1e-3).unwrap();
/// assert_eq!(const_step.current_time(), 0f64);
/// assert_eq!(const_step.next(), Some(0f64));
/// assert_eq!(const_step.next(), Some(1e-3));
/// assert_eq!(const_step.next(), Some(2e-3));
/// ```
///
/// ```ignore
/// # use moldybrody::approx::time::ConstStep;
/// # use moldybrody::approx::TimeIterator;
/// let mut const_step = ConstStep::new(1e-3).unwrap();
/// for time in const_step{
///     // something in here with time
/// }
/// ```
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ConstStep<T> {
    pub dt: T,
    /// step size of time iteration
    pub tmax: T,
    /// Maximum time
    current: T, //  current time of iterator
}

macro_rules! impl_const_step {
    ($ty : ident) => {
        impl ConstStep<$ty> {
            /// Return a new ConstStep with given time step
            ///
            /// # Examples
            /// ```
            /// # use moldybrody::approx::time::ConstStep;
            /// # use moldybrody::approx::TimeIterator;
            /// # use assert_approx_eq::assert_approx_eq;
            /// let const_step = ConstStep::new(1e-3).unwrap();
            /// assert_approx_eq!(const_step.dt, 1e-3);
            /// assert_approx_eq!(const_step.current, 0f64);
            /// assert_approx_eq!(const_step.tmax, std::f64::MAX);
            /// ```
            #[allow(dead_code)]
            pub fn new(dt: $ty) -> Result<Self, Error> {
                if dt < 1e-15 as $ty {
                    return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
                }
                Ok(Self {
                    current: 0 as $ty,
                    dt: dt,
                    tmax: std::$ty::MAX,
                })
            }

            #[allow(dead_code)]
            pub fn renew(&mut self) {
                self.current = 0 as $ty;
            }
        }

        impl TimeIterator<$ty> for ConstStep<$ty> {
            /// # Examples
            /// ```
            /// # use moldybrody::approx::time::ConstStep;
            /// # use moldybrody::approx::TimeIterator;
            /// # use assert_approx_eq::assert_approx_eq;
            /// let const_step = ConstStep::new(1e-3).unwrap();
            /// assert_approx_eq!(const_step.current_time(), 0f64);
            /// ```
            fn current_time(&self) -> $ty {
                self.current
            }

            /// # Examples
            /// ```
            /// # use moldybrody::approx::time::ConstStep;
            /// # use moldybrody::approx::TimeIterator;
            /// # use assert_approx_eq::assert_approx_eq;
            /// let const_step = ConstStep::new(1e-3).unwrap();
            /// assert_approx_eq!(const_step.dt(), 1e-3);
            /// ```
            fn dt(&self) -> $ty {
                self.dt
            }

            /// # Examples
            /// ```
            /// # use moldybrody::approx::time::ConstStep;
            /// # use moldybrody::approx::TimeIterator;
            /// # use assert_approx_eq::assert_approx_eq;
            /// let mut const_step = ConstStep::new(1e-3).unwrap();
            /// assert_approx_eq!(const_step.tmax, std::f64::MAX);
            ///
            /// const_step.set_tmax(1f64);
            /// assert_approx_eq!(const_step.tmax, 1f64);
            ///
            /// let mut max_t = 0f64;
            /// for t in const_step{
            ///     if t > max_t{
            ///         max_t = t;
            ///     }
            /// }
            /// max_t += 1e-3;
            /// assert_approx_eq!(max_t, 1f64);
            /// ```
            fn set_tmax(&mut self, tmax: $ty) -> Result<(), Error> {
                if tmax < 0 as $ty {
                    return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
                } else if tmax < 1e-10 as $ty {
                    self.tmax = std::$ty::MAX;
                } else {
                    self.tmax = tmax + self.dt * (0.5 as $ty);
                }

                Ok(())
            }

            fn into_diff(&self) -> TimeDiffIterator<Self> {
                TimeDiffIterator {
                    timeiter: self.clone(),
                }
            }
        }

        impl Iterator for ConstStep<$ty> {
            type Item = $ty;

            fn next(&mut self) -> Option<Self::Item> {
                let time = self.current;
                if time <= self.tmax {
                    self.current += self.dt;
                    return Some(time);
                } else {
                    return None;
                }
            }
        }
    };
}

impl_const_step!(f32);
impl_const_step!(f64);

// =============================================================================

/// Time iterator of which time step size is increasing exponentially
///
/// 초반의 time step을 매우 작은 단위에서부터 지수함수적으로 증가시키는 방식의 time iterator
///
/// # Examples
/// ```
/// # use moldybrody::approx::time::ExponentialStep;
/// # use moldybrody::approx::TimeIterator;
/// let mut exp_step = ExponentialStep::new(1e-5, 1e-3, 1).unwrap();
/// exp_step.set_inc(10f64).unwrap();
/// assert_eq!(exp_step.current_time(), 0f64);
/// assert_eq!(exp_step.next(), Some(0f64));
/// assert_eq!(exp_step.next(), Some(1e-5));
/// assert_eq!(exp_step.next(), Some(1.1e-4));
/// assert_eq!(exp_step.next(), Some(1.11e-3));
/// assert_eq!(exp_step.next(), Some(2.11e-3));
/// ```
///
/// ```ignore
/// # use moldybrody::approx::time::ExponentialStep;
/// # use moldybrody::approx::TimeIterator;
/// let mut exp_step = ExponentialStep::new(1e-5, 1e-3, 10).unwrap();
/// for time in exp_step{
///     // something in here with time
/// }
/// ```
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ExponentialStep<T> {
    pub dt_min: T,
    /// initial time step size
    pub dt_max: T,
    /// final time step size
    pub tmax: T,
    /// maximum time for simulation
    pub length: usize,
    /// increase period of step size
    current: T, //  current time
    dt: T,        //  current time step
    inc: T,       //  increase ratio of step size
    count: usize, //  count for next time step
}

macro_rules! impl_exp_step {
    ($ty : ident) => {
        impl ExponentialStep<$ty> {
            #[allow(dead_code)]
            pub fn new(dt_min: $ty, dt_max: $ty, length: usize) -> Result<Self, Error> {
                if dt_min < 1e-10 as $ty || dt_min > dt_max || length == 0 {
                    return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
                }

                let inc: $ty;
                if length == 1 {
                    inc = 1.2 as $ty;
                } else if length < 20 {
                    inc = (length as $ty) / (2 as $ty);
                } else {
                    inc = 10 as $ty;
                }

                Ok(Self {
                    dt_min: dt_min,
                    dt_max: dt_max,
                    tmax: std::$ty::MAX,
                    length: length,
                    current: 0 as $ty,
                    dt: dt_min,
                    inc: inc,
                    count: 0,
                })
            }

            #[allow(dead_code)]
            pub fn renew(&mut self) {
                self.current = 0 as $ty;
                self.dt = self.dt_min;
                self.count = 0;
            }

            #[allow(dead_code)]
            pub fn set_inc(&mut self, inc: $ty) -> Result<(), Error> {
                if inc <= 1 as $ty {
                    return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
                }
                self.inc = inc;
                Ok(())
            }
        }

        impl TimeIterator<$ty> for ExponentialStep<$ty> {
            fn current_time(&self) -> $ty {
                self.current
            }

            fn dt(&self) -> $ty {
                self.dt
            }

            fn set_tmax(&mut self, tmax: $ty) -> Result<(), Error> {
                if tmax < 0 as $ty {
                    return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
                } else if tmax < 1e-15 as $ty {
                    self.tmax = std::$ty::MAX;
                } else {
                    self.tmax = tmax + self.dt_max * (0.5 as $ty);
                }
                Ok(())
            }

            fn into_diff(&self) -> TimeDiffIterator<Self> {
                TimeDiffIterator {
                    timeiter: self.clone(),
                }
            }
        }

        impl Iterator for ExponentialStep<$ty> {
            type Item = $ty;

            fn next(&mut self) -> Option<Self::Item> {
                if self.dt < self.dt_max && self.count == self.length {
                    self.count = 0;
                    self.dt *= self.inc;
                    if self.dt > self.dt_max {
                        self.dt = self.dt_max;
                    }
                }
                let time = self.current;
                if time <= self.tmax {
                    self.current += self.dt;
                    self.count += 1;
                    return Some(time);
                } else {
                    return None;
                }
            }
        }
    };
}

impl_exp_step!(f32);
impl_exp_step!(f64);

#[cfg(test)]
mod test {
    use super::*;
    use serde_json::{from_str, to_string};
    #[test]
    fn test_serde_timeiter() {
        // Constant Step
        let mut cs = ConstStep::<f64>::new(1e-3).unwrap();
        cs.set_tmax(1.0).unwrap();
        let expected = r#"{"dt":0.001,"tmax":1.0005,"current":0.0}"#;
        assert_eq!(expected, to_string(&cs).unwrap());

        let expected: ConstStep<f64> = from_str(&expected).unwrap();
        assert_eq!(cs, expected);

        // Exponential Step
        let mut exps = ExponentialStep::<f64>::new(1e-3, 1e-2, 10).unwrap();
        exps.set_tmax(1.0).unwrap();
        let expected = r#"{"dt_min":0.001,"dt_max":0.01,"tmax":1.005,"length":10,"current":0.0,"dt":0.001,"inc":5.0,"count":0}"#;
        assert_eq!(expected, to_string(&exps).unwrap());

        let expected: ExponentialStep<f64> = from_str(&expected).unwrap();
        assert_eq!(exps, expected);
    }
}
