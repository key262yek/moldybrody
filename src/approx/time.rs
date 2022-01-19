//! Implement widely used time iterators
//!
//! 많은 영역에서 사용되고 있는 time iterator들을 구현해놓았습니다.
//! - Constant Step : 일정한 time step으로 시뮬레이션할 때 사용할 수 있습니다.
//! - Exponential Step : 시스템 초반에는 짧은 time step으로 시뮬레이션하고, 그 이후부터는 constant step으로 시뮬레이션 할 때 사용할 수 있습니다.
//! - Nearest Time : Run and Tumble 방식으로 움직이는 입자들의 경우, Run 과정에서는 굳이 짧은 시간 단위로 끊지 않아도 오차 없이 입자들의 운동을 기술할 수 있습니다. 따라서 우리는 각 입자의 Tumble time을 ordered list에 보관하고 가장 근접한 시점으로 바로 이동하는 방식을 채택할 수 있습니다.

use crate::prelude::*;

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
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct ConstStep{
    pub dt      : f64,          /// step size of time iteration
    pub tmax    : f64,          /// Maximum time
    current     : f64,          //  current time of iterator
}

impl ConstStep{
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
    pub fn new(dt : f64) -> Result<Self>{
        if dt < 1e-15{
            return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
        }
        Ok(Self{
            current : 0f64,
            dt      : dt,
            tmax    : std::f64::MAX,
        })
    }
}

impl TimeIterator<f64> for ConstStep{
    /// # Examples
    /// ```
    /// # use moldybrody::approx::time::ConstStep;
    /// # use moldybrody::approx::TimeIterator;
    /// # use assert_approx_eq::assert_approx_eq;
    /// let const_step = ConstStep::new(1e-3).unwrap();
    /// assert_approx_eq!(const_step.current_time(), 0f64);
    /// ```
    fn current_time(&self) -> f64{
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
    fn dt(&self) -> f64{
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
    fn set_tmax(&mut self, tmax : f64) -> Result<()>{
        if tmax < 0f64{
            return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
        }
        else if tmax < 1e-10{
            self.tmax = std::f64::MAX;
        }
        else{
            self.tmax = tmax + self.dt * 0.5;
        }

        Ok(())
    }

    fn into_diff(&self) -> TimeDiffIterator<Self>{
        TimeDiffIterator{
            timeiter : self.clone(),
        }
    }
}

impl Iterator for ConstStep{
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item>{
        let time = self.current;
        if time <= self.tmax{
            self.current += self.dt;
            return Some(time);
        }
        else{
            return None;
        }
    }
}


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
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct ExponentialStep{
    pub dt_min  : f64,                  /// initial time step size
    pub dt_max  : f64,                  /// final time step size
    pub tmax    : f64,                  /// maximum time for simulation
    pub length  : usize,                /// increase period of step size
    current : f64,                      //  current time
    dt      : f64,                      //  current time step
    inc     : f64,                      //  increase ratio of step size
    count   : usize,                    //  count for next time step
}

impl ExponentialStep{
    #[allow(dead_code)]
    pub fn new(dt_min : f64, dt_max : f64, length : usize) -> Result<Self>{
        if dt_min < 1e-10 || dt_min > dt_max || length == 0{
            return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
        }

        let inc : f64;
        if length == 1{
            inc = 1.2f64;
        }
        else if length < 20{
            inc = (length as f64) / 2f64;
        }
        else{
            inc = 10f64;       
        }

        Ok(Self{
            dt_min  : dt_min,
            dt_max  : dt_max,
            tmax    : std::f64::MAX,
            length  : length,
            current : 0f64,
            dt      : dt_min,
            inc     : inc,
            count   : 0,
        })
    }

    #[allow(dead_code)]
    pub fn set_inc(&mut self, inc : f64) -> Result<()>{
        if inc <= 1f64{
            return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
        }
        self.inc = inc;
        Ok(())
    }
}


impl TimeIterator<f64> for ExponentialStep{
    fn current_time(&self) -> f64{
        self.current
    }

    fn dt(&self) -> f64{
        self.dt
    }

    fn set_tmax(&mut self, tmax : f64) -> Result<()>{
        if tmax < 0f64{
            return Err(Error::make_error_syntax(ErrorCode::InvalidArgumentInput));
        }
        else if tmax < 1e-15{
            self.tmax = std::f64::MAX;
        }
        else{
            self.tmax = tmax;
        }
        Ok(())
    }

    fn into_diff(&self) -> TimeDiffIterator<Self>{
        TimeDiffIterator{
            timeiter : self.clone(),
        }
    }
}

impl Iterator for ExponentialStep{
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item>{
        if self.dt < self.dt_max && self.count == self.length{
            self.count = 0;
            self.dt *= self.inc;
            if self.dt > self.dt_max{
                self.dt = self.dt_max;
            }
        }
        let time = self.current;
        if time <= self.tmax{
            self.current += self.dt;
            self.count += 1;
            return Some(time);
        }
        else{
            return None;
        }
    }
}


