
use crate::error::Error;
use approx::AbsDiffEq;

pub struct ButcherTableau<T>{
    aij : Vec<Vec<T>>,
    bi : Vec<T>,
    ci : Vec<T>,
}

pub(crate) enum NewtonButcherTableauBuilder<T>{
    Euler,
    MidPoint,
    Ssprk3,
    Classic,
    Heun(usize),
    Ralston(usize),
    Generic(usize, T),
}

impl<T> NewtonButcherTableauBuilder<T>{
    pub(crate) fn new(name : Option<&str>, order : Option<usize>, alpha : Option<T>) -> Result<Self, Error>{
        match name{
            Some("Euler") => Ok(NewtonButcherTableauBuilder::Euler),
            Some("MidPoint") => Ok(NewtonButcherTableauBuilder::MidPoint),
            Some("SSPRK3") => Ok(NewtonButcherTableauBuilder::Ssprk3),
            Some("Classic") => Ok(NewtonButcherTableauBuilder::Classic),
            Some("Heun") => {
                match order {
                    Some(v) => Ok(NewtonButcherTableauBuilder::Heun(v)),
                    None =>  Err(Error::make_error_msg("Heun method requires order attribute.".to_string())),
                }
            },
            Some("Ralston") => {
                match order {
                    Some(v) => Ok(NewtonButcherTableauBuilder::Ralston(v)),
                    None => Err(Error::make_error_msg("Ralston method requires order attribute.".to_string())),
                }
            },
            Some("Generic") => {
                match (order, alpha) {
                    (Some(v), Some(a)) => Ok(NewtonButcherTableauBuilder::Generic(v, a)),
                    (Some(_), None) => Err(Error::make_error_msg("Generic Runge-Kutta method requires parameter alpha in [0, 1]".to_string())),
                    (None, _) => Err(Error::make_error_msg("Generic Runge-Kutta method requires order".to_string())),
                }
            }
            None | Some(_) => Err(Error::make_error_msg("Butcher Tableau builder requires valid name of tableau. ex) Euler, Midpoint, SSPRK3, Classic, Heun, Ralston, Generic".to_string())),
        }
    }
}

macro_rules! impl_bt_builder{
    ($ty : ident) => {
        impl NewtonButcherTableauBuilder<$ty> {
            pub(crate) fn build(self) -> Result<ButcherTableau<$ty>, Error>{
                match self{
                    NewtonButcherTableauBuilder::Euler => {
                        Ok(ButcherTableau{
                            aij : vec![vec![0.]],
                            bi : vec![1.],
                            ci : vec![0.],
                        })
                    },
                    NewtonButcherTableauBuilder::MidPoint => {
                        Ok(ButcherTableau{
                            aij : vec![vec![0., 0.], vec![0.5, 0.]],
                            bi : vec![0., 1.],
                            ci : vec![0., 0.5],
                        })
                    },
                    NewtonButcherTableauBuilder::Ssprk3 => {
                        let t = 1. / 6.;
                        Ok(ButcherTableau{
                            aij : vec![vec![0., 0., 0.], vec![1., 0., 0.], vec![0.25, 0.25, 0.]],
                            bi : vec![t, t, 4. * t],
                            ci : vec![0., 1., 0.5],
                        })  
                    },
                    NewtonButcherTableauBuilder::Classic => {
                        let t = 1. / 6.;
                        Ok(ButcherTableau{
                            aij : vec![vec![0., 0., 0., 0.], vec![0.5, 0., 0., 0.], vec![0., 0.5, 0., 0.], vec![0., 0., 1., 0.]],
                            bi : vec![t, 2. * t, 2. * t, t],
                            ci : vec![0., 0.5, 0.5, 1.],
                        })
                    },
                    NewtonButcherTableauBuilder::Heun(order) => {
                        match order {
                            2 => {
                                Ok(ButcherTableau{
                                    aij : vec![vec![0., 0.], vec![1., 0.]],
                                    bi : vec![0.5, 0.5],
                                    ci : vec![0., 1.],
                                })
                            },
                            3 => {
                                let t = 1. / 3.;
                                Ok(ButcherTableau{
                                    aij : vec![vec![0., 0., 0.], vec![t, 0., 0.], vec![0., 2. * t, 0.]],
                                    bi : vec![0.25, 0., 0.75],
                                    ci : vec![0., t, 2. * t],
                                })
                            }, 
                            _ => {
                                Err(Error::make_error_msg(("Heun methods exist only for 2nd and 3rd order".to_string())))
                            }
                        }
                    },
                    NewtonButcherTableauBuilder::Ralston(order) => {
                        match order {
                            2 => {
                                Ok(ButcherTableau{
                                    aij : vec![vec![0., 0.], vec![2. / 3., 0.]],
                                    bi : vec![0.25, 0.75],
                                    ci : vec![0., 2. / 3.],
                                })
                            },
                            3 => {
                                let t = 1. / 9.;
                                Ok(ButcherTableau{
                                    aij : vec![vec![0., 0., 0.], vec![0.5, 0., 0.], vec![0., 0.75, 0.]],
                                    bi : vec![2. * t, 3. * t, 4. * t],
                                    ci : vec![0., 0.5, 0.75],
                                })  
                            }, 
                            _ => {
                                Err(Error::make_error_msg(("Raltson methods exist only for 2nd and 3rd order".to_string())))
                            }
                        }
                    },
                    NewtonButcherTableauBuilder::Generic(order, alpha) => {
                        match order {
                            2 => {
                                let t = 1. / (2. * alpha);
                                Ok(ButcherTableau{
                                    aij : vec![vec![0., 0.], vec![alpha, 0.]],
                                    bi : vec![1. - t, t],
                                    ci : vec![0., alpha],
                                })
                            },
                            3 => {
                                if alpha.abs_diff_eq(&0., 1e-5) || alpha.abs_diff_eq(&(2. / 3.), 1e-5) || alpha.abs_diff_eq(&1., 1e-5){
                                    return Err(Error::make_error_msg(("Generic 3rd order Runge-Kutta method cannot be made with alpha = 0, 2/3, 1".to_string())))
                                }
                                let t = (1. - alpha) / (3. * alpha - 2.);
                                Ok(ButcherTableau{
                                    aij : vec![vec![0., 0., 0.], vec![alpha, 0., 0.], vec![1. + t / alpha, - t / alpha, 0.]],
                                    bi : vec![0.5 - 1. / (6. * alpha), 1. / (6. * alpha * (1. - alpha)), - t / 6.],
                                    ci : vec![0., alpha, 1.],
                                })
                            }, 
                            _ => {
                                Err(Error::make_error_msg(("Generic Runge-Kutta tableau exist only for 2nd and 3rd order".to_string())))
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


#[cfg(test)]
mod test{
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_tableau_builder(){
        let euler : ButcherTableau<f32> = NewtonButcherTableauBuilder::<f32>::new(Some("Euler"), None, None).unwrap()
                                            .build().unwrap();
        assert_abs_diff_eq!(euler.aij[0][0], 0.0);
        assert_abs_diff_eq!(euler.bi[0], 1.);
        assert_abs_diff_eq!(euler.ci[0], 0.);
    }
}