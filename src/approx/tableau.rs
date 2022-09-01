
use crate::error::Error;
use crate::argument::CommandBuilder;
use approx::AbsDiffEq;
use clap::Arg;
use clap::ArgMatches;


pub struct ButcherTableau<T>{
    aij : Vec<Vec<T>>,
    bi : Vec<T>,
    ci : Vec<T>,
}

pub trait ApproxMethod<T>{
    fn tableau(&self) -> &ButcherTableau<T>;
}

pub struct NewtonEulerMethod<T>{
    tableau : ButcherTableau<T>
}

impl<T> ApproxMethod<T> for NewtonEulerMethod<T>{
    fn tableau(&self) -> &ButcherTableau<T>{
        &self.tableau
    }
}

macro_rules! impl_euler{
    ($ty : ident) => {
        impl NewtonEulerMethod<$ty>{
            pub fn new() -> Self{
                Self {
                    tableau : ButcherTableau{
                        aij : vec![vec![0.]],
                        bi : vec![1.],
                        ci : vec![0.],
                    }
                }
            }
        }
        
        impl<'h> CommandBuilder<'h, 0> for NewtonEulerMethod<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }
        
        impl From<&ArgMatches> for NewtonEulerMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                NewtonEulerMethod::<$ty>::new()
            }
        }
    }
}

impl_euler!(f32);
impl_euler!(f64);

pub struct NewtonMidPointMethod<T>{
    tableau : ButcherTableau<T>
}

impl<T> ApproxMethod<T> for NewtonMidPointMethod<T>{
    fn tableau(&self) -> &ButcherTableau<T>{
        &self.tableau
    }
}

macro_rules! impl_midpoint{
    ($ty : ident) => {
        impl NewtonMidPointMethod<$ty>{
            pub fn new() -> Self{
                Self {
                    tableau : ButcherTableau{
                        aij : vec![vec![0., 0.], vec![0.5, 0.]],
                        bi : vec![0., 1.],
                        ci : vec![0., 0.5],
                    }
                }
            }
        }
        
        impl<'h> CommandBuilder<'h, 0> for NewtonMidPointMethod<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }
        
        impl From<&ArgMatches> for NewtonMidPointMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                NewtonMidPointMethod::<$ty>::new()
            }
        }
    }
}

impl_midpoint!(f32);
impl_midpoint!(f64);

pub struct NewtonSSPRK3Method<T>{
    tableau : ButcherTableau<T>
}

impl<T> ApproxMethod<T> for NewtonSSPRK3Method<T>{
    fn tableau(&self) -> &ButcherTableau<T>{
        &self.tableau
    }
}

macro_rules! impl_ssprk3{
    ($ty : ident) => {
        impl NewtonSSPRK3Method<$ty>{
            pub fn new() -> Self{
                let t = 1. / 6.;
                Self {
                    tableau : ButcherTableau{
                        aij : vec![vec![0., 0., 0.], vec![1., 0., 0.], vec![0.25, 0.25, 0.]],
                        bi : vec![t, t, 4. * t],
                        ci : vec![0., 1., 0.5],
                    }
                }
            }
        }
        
        impl<'h> CommandBuilder<'h, 0> for NewtonSSPRK3Method<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }
        
        impl From<&ArgMatches> for NewtonSSPRK3Method<$ty> {
            fn from(m: &ArgMatches) -> Self {
                NewtonSSPRK3Method::<$ty>::new()
            }
        }
    }
}

impl_ssprk3!(f32);
impl_ssprk3!(f64);


pub struct NewtonClassicMethod<T>{
    tableau : ButcherTableau<T>
}

impl<T> ApproxMethod<T> for NewtonClassicMethod<T>{
    fn tableau(&self) -> &ButcherTableau<T>{
        &self.tableau
    }
}

macro_rules! impl_classic{
    ($ty : ident) => {
        impl NewtonClassicMethod<$ty>{
            pub fn new() -> Self{
                let t = 1. / 6.;
                Self {
                    tableau : ButcherTableau{
                        aij : vec![vec![0., 0., 0., 0.], vec![0.5, 0., 0., 0.], vec![0., 0.5, 0., 0.], vec![0., 0., 1., 0.]],
                        bi : vec![t, 2. * t, 2. * t, t],
                        ci : vec![0., 0.5, 0.5, 1.],
                    }
                }
            }
        }
        
        impl<'h> CommandBuilder<'h, 0> for NewtonClassicMethod<$ty> {
            fn args() -> [Arg<'h>; 0] {
                []
            }
        }
        
        impl From<&ArgMatches> for NewtonClassicMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                NewtonClassicMethod::<$ty>::new()
            }
        }
    }
}

impl_classic!(f32);
impl_classic!(f64);

pub struct NewtonHeunMethod<T>{
    tableau : ButcherTableau<T>
}

impl<T> ApproxMethod<T> for NewtonHeunMethod<T>{
    fn tableau(&self) -> &ButcherTableau<T>{
        &self.tableau
    }
}

macro_rules! impl_heun{
    ($ty : ident) => {
        impl NewtonHeunMethod<$ty>{
            pub fn new(order : usize) -> Self{
                let tableau = match order {
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
                        panic!("Heun methods exist only for 2nd and 3rd order")
                    }
                };
                Self{ tableau }
            }
        }
        
        impl<'h> CommandBuilder<'h, 1> for NewtonHeunMethod<$ty> {
            fn args() -> [Arg<'h>; 1] {
                [
                    Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Heun method. You can choose between 2 or 3. Default value = 3")
                ]
            }
        }
        
        impl From<&ArgMatches> for NewtonHeunMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                let order : usize = *m.get_one::<usize>("order").unwrap();
                NewtonHeunMethod::<$ty>::new(order)
            }
        }
    }
}

impl_heun!(f32);
impl_heun!(f64);

pub struct NewtonRalstonMethod<T>{
    tableau : ButcherTableau<T>
}

impl<T> ApproxMethod<T> for NewtonRalstonMethod<T>{
    fn tableau(&self) -> &ButcherTableau<T>{
        &self.tableau
    }
}

macro_rules! impl_raltson{
    ($ty : ident) => {
        impl NewtonRalstonMethod<$ty>{
            pub fn new(order : usize) -> Self{
                let tableau = match order {
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
                        panic!("Raltson methods exist only for 2nd and 3rd order")
                    }
                };
                Self{ tableau }
            }
        }
        
        impl<'h> CommandBuilder<'h, 1> for NewtonRalstonMethod<$ty> {
            fn args() -> [Arg<'h>; 1] {
                [
                    Arg::new("order")
                    .long("order")
                    .value_name("ORDER")
                    .default_value("3")
                    .value_parser(clap::value_parser!(usize))
                    .help("Order of Raltson method. You can choose between 2 or 3. Default value = 3")
                ]
            }
        }
        
        impl From<&ArgMatches> for NewtonRalstonMethod<$ty> {
            fn from(m: &ArgMatches) -> Self {
                let order : usize = *m.get_one::<usize>("order").unwrap();
                NewtonRalstonMethod::<$ty>::new(order)
            }
        }
    }
}

impl_raltson!(f32);
impl_raltson!(f64);

pub struct NewtonGenericRungeKuttaMethod<T>{
    tableau : ButcherTableau<T>
}

impl<T> ApproxMethod<T> for NewtonGenericRungeKuttaMethod<T>{
    fn tableau(&self) -> &ButcherTableau<T>{
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
                            bi : vec![1. - t, t],
                            ci : vec![0., alpha],
                        }
                    },
                    3 => {
                        if alpha.abs_diff_eq(&0., 1e-5) || alpha.abs_diff_eq(&(2. / 3.), 1e-5) || alpha.abs_diff_eq(&1., 1e-5){
                            return panic!("Generic 3rd order Runge-Kutta method cannot be made with alpha = 0, 2/3, 1")
                        }
                        let t = (1. - alpha) / (3. * alpha - 2.);
                        ButcherTableau{
                            aij : vec![vec![0., 0., 0.], vec![alpha, 0., 0.], vec![1. + t / alpha, - t / alpha, 0.]],
                            bi : vec![0.5 - 1. / (6. * alpha), 1. / (6. * alpha * (1. - alpha)), - t / 6.],
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
mod test{
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_tableau_builder(){
        // let euler : ButcherTableau<f32> = NewtonButcherTableauBuilder::<f32>::new(Some("Euler"), None, None).unwrap()
        //                                     .build().unwrap();
        // assert_abs_diff_eq!(euler.aij[0][0], 0.0);
        // assert_abs_diff_eq!(euler.bi[0], 1.);
        // assert_abs_diff_eq!(euler.ci[0], 0.);
    }
}