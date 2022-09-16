use crate::argument::CommandBuilder;
use approx::AbsDiffEq;
use clap::Arg;
use clap::ArgMatches;
use std::fmt::Debug;

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct ButcherTableau<T> {
    pub(crate) aij: Vec<Vec<T>>,
    pub(crate) bj: Vec<T>,
    pub(crate) ci: Vec<T>,
}

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct TableauIterator<'a, T> {
    pub(crate) tableau : &'a ButcherTableau<T>,
    pub(crate) idx : usize
}

impl<T> ButcherTableau<T> {
    pub fn len(&self) -> usize {
        self.bj.len()
    }

    pub fn into_iter(&self) -> TableauIterator<T>{
        TableauIterator{
            tableau : &self,
            idx : 0,
        }
    }
}


impl<T> AbsDiffEq for ButcherTableau<T>
where
    T: AbsDiffEq,
    <T as AbsDiffEq>::Epsilon: Copy,
{
    type Epsilon = <T as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        <T as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let (n, m) = (self.ci.len(), self.bj.len());
        let (n2, m2) = (other.ci.len(), other.bj.len());

        if n != n2 || m != m2 {
            return false;
        }

        for i in 0..n {
            if !self.ci[i].abs_diff_eq(&other.ci[i], epsilon) {
                return false;
            }

            for j in 0..m {
                if !self.aij[i][j].abs_diff_eq(&other.aij[i][j], epsilon) {
                    return false;
                }
            }
        }

        for j in 0..m {
            if !self.bj[j].abs_diff_eq(&other.bj[j], epsilon) {
                return false;
            }
        }

        return true;
    }
}

impl<'a, T> Iterator for TableauIterator<'a, T>
    where T : num_traits::One + Copy{
    type Item = (T, &'a Vec<T>);

    fn next(&mut self) -> Option<Self::Item>{
        let length = self.tableau.len();
        if self.idx < length{
            let temp = self.idx;
            self.idx += 1;
            return Some((self.tableau.ci[temp], &self.tableau.aij[temp]));
        } else if self.idx == length {
            self.idx += 1;
            return Some((<T as num_traits::One>::one(), &self.tableau.bj));
        } else {
            return None;
        }
    }
}

pub trait ApproxMethod<T> : Debug {
    fn tableau(&self) -> &ButcherTableau<T>;
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NewtonEulerMethod<T> {
    tableau: ButcherTableau<T>,
}

impl<T> ApproxMethod<T> for NewtonEulerMethod<T> 
    where T : Debug{
    fn tableau(&self) -> &ButcherTableau<T> {
        &self.tableau
    }
}

macro_rules! impl_euler {
    ($ty : ident) => {
        impl NewtonEulerMethod<$ty> {
            pub fn new() -> Self {
                Self {
                    tableau: ButcherTableau {
                        aij: vec![vec![0.]],
                        bj: vec![1.],
                        ci: vec![0.],
                    },
                }
            }
        }

        impl<'h> CommandBuilder<'h, 0> for NewtonEulerMethod<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }

        impl From<&ArgMatches> for NewtonEulerMethod<$ty> {
            fn from(_m: &ArgMatches) -> Self {
                NewtonEulerMethod::<$ty>::new()
            }
        }
    };
}

impl_euler!(f32);
impl_euler!(f64);

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NewtonMidPointMethod<T> {
    tableau: ButcherTableau<T>,
}

impl<T> ApproxMethod<T> for NewtonMidPointMethod<T> 
    where T : Debug{
    fn tableau(&self) -> &ButcherTableau<T> {
        &self.tableau
    }
}

macro_rules! impl_midpoint {
    ($ty : ident) => {
        impl NewtonMidPointMethod<$ty> {
            pub fn new() -> Self {
                Self {
                    tableau: ButcherTableau {
                        aij: vec![vec![0., 0.], vec![0.5, 0.]],
                        bj: vec![0., 1.],
                        ci: vec![0., 0.5],
                    },
                }
            }
        }

        impl<'h> CommandBuilder<'h, 0> for NewtonMidPointMethod<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }

        impl From<&ArgMatches> for NewtonMidPointMethod<$ty> {
            fn from(_m: &ArgMatches) -> Self {
                NewtonMidPointMethod::<$ty>::new()
            }
        }
    };
}

impl_midpoint!(f32);
impl_midpoint!(f64);

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NewtonSSPRK3Method<T> {
    tableau: ButcherTableau<T>,
}

impl<T> ApproxMethod<T> for NewtonSSPRK3Method<T> 
    where T : Debug{
    fn tableau(&self) -> &ButcherTableau<T> {
        &self.tableau
    }
}

macro_rules! impl_ssprk3 {
    ($ty : ident) => {
        impl NewtonSSPRK3Method<$ty> {
            pub fn new() -> Self {
                let t = 1. / 6.;
                Self {
                    tableau: ButcherTableau {
                        aij: vec![vec![0., 0., 0.], vec![1., 0., 0.], vec![0.25, 0.25, 0.]],
                        bj: vec![t, t, 4. * t],
                        ci: vec![0., 1., 0.5],
                    },
                }
            }
        }

        impl<'h> CommandBuilder<'h, 0> for NewtonSSPRK3Method<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }

        impl From<&ArgMatches> for NewtonSSPRK3Method<$ty> {
            fn from(_m: &ArgMatches) -> Self {
                NewtonSSPRK3Method::<$ty>::new()
            }
        }
    };
}

impl_ssprk3!(f32);
impl_ssprk3!(f64);

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NewtonClassicMethod<T> {
    tableau: ButcherTableau<T>,
}

impl<T> ApproxMethod<T> for NewtonClassicMethod<T> 
    where T : Debug{
    fn tableau(&self) -> &ButcherTableau<T> {
        &self.tableau
    }
}

macro_rules! impl_classic {
    ($ty : ident) => {
        impl NewtonClassicMethod<$ty> {
            pub fn new() -> Self {
                let t = 1. / 6.;
                Self {
                    tableau: ButcherTableau {
                        aij: vec![
                            vec![0., 0., 0., 0.],
                            vec![0.5, 0., 0., 0.],
                            vec![0., 0.5, 0., 0.],
                            vec![0., 0., 1., 0.],
                        ],
                        bj: vec![t, 2. * t, 2. * t, t],
                        ci: vec![0., 0.5, 0.5, 1.],
                    },
                }
            }
        }

        impl<'h> CommandBuilder<'h, 0> for NewtonClassicMethod<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }

        impl From<&ArgMatches> for NewtonClassicMethod<$ty> {
            fn from(_m: &ArgMatches) -> Self {
                NewtonClassicMethod::<$ty>::new()
            }
        }
    };
}

impl_classic!(f32);
impl_classic!(f64);

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NewtonHeunMethod<T> {
    tableau: ButcherTableau<T>,
}

impl<T> ApproxMethod<T> for NewtonHeunMethod<T> 
    where T : Debug{
    fn tableau(&self) -> &ButcherTableau<T> {
        &self.tableau
    }
}

macro_rules! impl_heun {
    ($ty : ident) => {
        impl NewtonHeunMethod<$ty> {
            pub fn new(order: usize) -> Self {
                let tableau = match order {
                    2 => ButcherTableau {
                        aij: vec![vec![0., 0.], vec![1., 0.]],
                        bj: vec![0.5, 0.5],
                        ci: vec![0., 1.],
                    },
                    3 => {
                        let t = 1. / 3.;
                        ButcherTableau {
                            aij: vec![vec![0., 0., 0.], vec![t, 0., 0.], vec![0., 2. * t, 0.]],
                            bj: vec![0.25, 0., 0.75],
                            ci: vec![0., t, 2. * t],
                        }
                    }
                    _ => {
                        panic!("Heun methods exist only for 2nd and 3rd order")
                    }
                };
                Self { tableau }
            }
        }

        impl<'h> CommandBuilder<'h, 1> for NewtonHeunMethod<$ty> {
            fn args() -> [Arg<'h>; 1] {
                [Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Heun method. You can choose between 2 or 3. Default value = 3")]
            }
        }

        impl From<&ArgMatches> for NewtonHeunMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                let order: usize = *m.get_one::<usize>("order").unwrap();
                NewtonHeunMethod::<$ty>::new(order)
            }
        }
    };
}

impl_heun!(f32);
impl_heun!(f64);

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NewtonRalstonMethod<T> {
    tableau: ButcherTableau<T>,
}

impl<T> ApproxMethod<T> for NewtonRalstonMethod<T> 
    where T : Debug{
    fn tableau(&self) -> &ButcherTableau<T> {
        &self.tableau
    }
}

macro_rules! impl_raltson {
    ($ty : ident) => {
        impl NewtonRalstonMethod<$ty> {
            pub fn new(order: usize) -> Self {
                let tableau = match order {
                    2 => ButcherTableau {
                        aij: vec![vec![0., 0.], vec![2. / 3., 0.]],
                        bj: vec![0.25, 0.75],
                        ci: vec![0., 2. / 3.],
                    },
                    3 => {
                        let t = 1. / 9.;
                        ButcherTableau {
                            aij: vec![vec![0., 0., 0.], vec![0.5, 0., 0.], vec![0., 0.75, 0.]],
                            bj: vec![2. * t, 3. * t, 4. * t],
                            ci: vec![0., 0.5, 0.75],
                        }
                    }
                    _ => {
                        panic!("Raltson methods exist only for 2nd and 3rd order")
                    }
                };
                Self { tableau }
            }
        }

        impl<'h> CommandBuilder<'h, 1> for NewtonRalstonMethod<$ty> {
            fn args() -> [Arg<'h>; 1] {
                [Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help(
                        "Order of Raltson method. You can choose between 2 or 3. Default value = 3",
                    )]
            }
        }

        impl From<&ArgMatches> for NewtonRalstonMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                let order: usize = *m.get_one::<usize>("order").unwrap();
                NewtonRalstonMethod::<$ty>::new(order)
            }
        }
    };
}

impl_raltson!(f32);
impl_raltson!(f64);

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NewtonGenericRungeKuttaMethod<T> {
    tableau: ButcherTableau<T>,
}

impl<T> ApproxMethod<T> for NewtonGenericRungeKuttaMethod<T> 
    where T : Debug{
    fn tableau(&self) -> &ButcherTableau<T> {
        &self.tableau
    }
}

macro_rules! impl_generic_runge_kutta{
    ($ty : ident) => {
        impl NewtonGenericRungeKuttaMethod<$ty>{
            pub fn new(order : usize, alpha : $ty) -> Self{
                let tableau = match order {
                    2 => {
                        let t = 1. / (2. * alpha);
                        ButcherTableau{
                            aij : vec![vec![0., 0.], vec![alpha, 0.]],
                            bj : vec![1. - t, t],
                            ci : vec![0., alpha],
                        }
                    },
                    3 => {
                        if alpha.abs_diff_eq(&0., 1e-5) || alpha.abs_diff_eq(&(2. / 3.), 1e-5) || alpha.abs_diff_eq(&1., 1e-5){
                            panic!("Generic 3rd order Runge-Kutta method cannot be made with alpha = 0, 2/3, 1")
                        }
                        let t = (1. - alpha) / (3. * alpha - 2.);
                        ButcherTableau{
                            aij : vec![vec![0., 0., 0.], vec![alpha, 0., 0.], vec![1. + t / alpha, - t / alpha, 0.]],
                            bj : vec![0.5 - 1. / (6. * alpha), 1. / (6. * alpha * (1. - alpha)), - t / 6.],
                            ci : vec![0., alpha, 1.],
                        }
                    },
                    _ => {
                        panic!("Raltson methods exist only for 2nd and 3rd order")
                    }
                };
                Self{ tableau }
            }
        }

        impl<'h> CommandBuilder<'h, 2> for NewtonGenericRungeKuttaMethod<$ty> {
            fn args() -> [Arg<'h>; 2] {
                [
                    Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of generic Runge-Kutta method. You can choose between 2 or 3. Default value = 3"),
                    Arg::new("alpha")
                    .long("alpha")
                    .value_name("ALPHA")
                    .default_value("0.5")
                    .value_parser(clap::value_parser!($ty))
                    .help("variable alpha for generic Runge-Kutta method.")
                ]
            }
        }

        impl From<&ArgMatches> for NewtonGenericRungeKuttaMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                let order : usize = *m.get_one::<usize>("order").unwrap();
                let alpha : $ty = *m.get_one::<$ty>("alpha").unwrap();
                NewtonGenericRungeKuttaMethod::<$ty>::new(order, alpha)
            }
        }
    }
}

impl_generic_runge_kutta!(f32);
impl_generic_runge_kutta!(f64);

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_abs_diff_eq;
    use clap::Command;

    #[test]
    fn test_tableau() {
        let c = Command::new("test").get_matches_from(vec!["test"]);
        let euler = NewtonEulerMethod::<f32>::from(&c);

        assert_abs_diff_eq!(
            euler.tableau,
            &ButcherTableau {
                aij: vec![vec![0.]],
                bj: vec![1.],
                ci: vec![0.],
            },
            epsilon = 1e-5
        );

        let midpoint = NewtonMidPointMethod::<f32>::from(&c);
        assert_abs_diff_eq!(
            midpoint.tableau,
            &ButcherTableau {
                aij: vec![vec![0., 0.], vec![0.5, 0.]],
                bj: vec![0., 1.],
                ci: vec![0., 0.5],
            },
            epsilon = 1e-5
        );

        let ssprk3 = NewtonSSPRK3Method::<f32>::from(&c);
        let t: f32 = 1. / 6.;
        assert_abs_diff_eq!(
            ssprk3.tableau,
            &ButcherTableau {
                aij: vec![vec![0., 0., 0.], vec![1., 0., 0.], vec![0.25, 0.25, 0.]],
                bj: vec![t, t, 4. * t],
                ci: vec![0., 1., 0.5],
            },
            epsilon = 1e-5
        );

        let classic = NewtonClassicMethod::<f32>::from(&c);
        assert_abs_diff_eq!(
            classic.tableau,
            &ButcherTableau {
                aij: vec![
                    vec![0., 0., 0., 0.],
                    vec![0.5, 0., 0., 0.],
                    vec![0., 0.5, 0., 0.],
                    vec![0., 0., 1., 0.],
                ],
                bj: vec![t, 2. * t, 2. * t, t],
                ci: vec![0., 0.5, 0.5, 1.],
            },
            epsilon = 1e-5
        );

        let c = Command::new("test")
            .arg(
                Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Heun method. You can choose between 2 or 3. Default value = 3"),
            )
            .get_matches_from(vec!["test", "--order", "3"]);
        let heun = NewtonHeunMethod::<f32>::from(&c);
        let t: f32 = 1. / 3.;
        assert_abs_diff_eq!(
            heun.tableau,
            &ButcherTableau {
                aij: vec![vec![0., 0., 0.], vec![t, 0., 0.], vec![0., 2. * t, 0.]],
                bj: vec![0.25, 0., 0.75],
                ci: vec![0., t, 2. * t],
            },
            epsilon = 1e-5
        );

        let raltson = NewtonRalstonMethod::<f32>::from(&c);
        let t: f32 = 1. / 9.;
        assert_abs_diff_eq!(
            raltson.tableau,
            &ButcherTableau {
                aij: vec![vec![0., 0., 0.], vec![0.5, 0., 0.], vec![0., 0.75, 0.]],
                bj: vec![2. * t, 3. * t, 4. * t],
                ci: vec![0., 0.5, 0.75],
            },
            epsilon = 1e-5
        );

        let c = Command::new("test")
            .arg(
                Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Heun method. You can choose between 2 or 3. Default value = 3"),
            )
            .arg(
                Arg::new("alpha")
                    .long("alpha")
                    .value_name("ALPHA")
                    .default_value("0.5")
                    .value_parser(clap::value_parser!(f32))
                    .help("variable alpha for generic Runge-Kutta method."),
            )
            .get_matches_from(vec!["test", "--order", "3", "--alpha", "0.5"]);
        let runge_kutta = NewtonGenericRungeKuttaMethod::<f32>::from(&c);
        let alpha: f32 = 0.5;
        let t: f32 = (1. - alpha) / (3. * alpha - 2.);
        assert_abs_diff_eq!(
            runge_kutta.tableau,
            &ButcherTableau {
                aij: vec![
                    vec![0., 0., 0.],
                    vec![alpha, 0., 0.],
                    vec![1. + t / alpha, -t / alpha, 0.]
                ],
                bj: vec![
                    0.5 - 1. / (6. * alpha),
                    1. / (6. * alpha * (1. - alpha)),
                    -t / 6.
                ],
                ci: vec![0., alpha, 1.],
            },
            epsilon = 1e-5
        );
    }

    #[test]
    #[should_panic]
    fn test_tableau_exception() {
        let c = Command::new("test")
            .arg(
                Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Heun method. You can choose between 2 or 3. Default value = 3"),
            )
            .arg(
                Arg::new("alpha")
                    .long("alpha")
                    .value_name("ALPHA")
                    .default_value("0.5")
                    .value_parser(clap::value_parser!(f32))
                    .help("variable alpha for generic Runge-Kutta method."),
            )
            .get_matches_from(vec!["test", "--order", "3", "--alpha", "0.66666666"]);
        let _runge_kutta = NewtonGenericRungeKuttaMethod::<f32>::from(&c);
    }

    #[test]
    #[should_panic]
    fn test_tableau_exception2() {
        let c = Command::new("test")
            .arg(
                Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Heun method. You can choose between 2 or 3. Default value = 3"),
            )
            .arg(
                Arg::new("alpha")
                    .long("alpha")
                    .value_name("ALPHA")
                    .default_value("0.5")
                    .value_parser(clap::value_parser!(f32))
                    .help("variable alpha for generic Runge-Kutta method."),
            )
            .get_matches_from(vec!["test", "--order", "3", "--alpha", "0.0"]);
        let _runge_kutta = NewtonGenericRungeKuttaMethod::<f32>::from(&c);
    }

    #[test]
    #[should_panic]
    fn test_tableau_exception3() {
        let c = Command::new("test")
            .arg(
                Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Heun method. You can choose between 2 or 3. Default value = 3"),
            )
            .arg(
                Arg::new("alpha")
                    .long("alpha")
                    .value_name("ALPHA")
                    .default_value("0.5")
                    .value_parser(clap::value_parser!(f32))
                    .help("variable alpha for generic Runge-Kutta method."),
            )
            .get_matches_from(vec!["test", "--order", "3", "--alpha", "1.0"]);
        let _runge_kutta = NewtonGenericRungeKuttaMethod::<f32>::from(&c);
    }
}
