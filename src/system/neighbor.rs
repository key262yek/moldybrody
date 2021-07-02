//! Specific implementation for structures with the Neighbor trait
//!
//! Neighbor trait을 가지는 여러 structure들을 정의하는 module입니다.

use crate::prelude::*;

/// Open-Ball neighborhood corresponds to [`Cartessian1D`](struct@Cartessian1D)
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct OpenBall1D<'a>{
    /// Point of ball center
    pub center : &'a Cartessian1D<f64>,
}

impl Neighbor<'_> for OpenBall1D<'_>{}

impl<'a> OpenBall1D<'a>{
    /// Initializing [`OpenBall1D`](struct@OpenBall1D) with a reference of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::Cartessian1D;
    /// # use moldybrody::system::neighbor::OpenBall1D;
    /// let v = Cartessian1D{coord : 0.0};
    /// let n = OpenBall1D::new(&v);
    /// ```
    pub fn new(center : &'a Cartessian1D<f64>) -> Self{
        OpenBall1D{
            center,
        }
    }
}

/// Open-Ball neighborhood corresponds to [`Cartessian1D`](struct@Cartessian1D)
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct OpenBallND<'a>{
    /// Point of ball center
    pub center : &'a CartessianND<f64>,
}

impl Neighbor<'_> for OpenBallND<'_>{}

impl<'a> OpenBallND<'a>{
    /// Initializing [`OpenBallND`](struct@OpenBallND) with a reference of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// # use moldybrody::system::point::CartessianND;
    /// # use moldybrody::system::neighbor::OpenBallND;
    /// let v = CartessianND{coord : vec![0.0; 4]};
    /// let n = OpenBallND::new(&v);
    /// ```
    pub fn new(center : &'a CartessianND<f64>) -> Self{
        OpenBallND{
            center,
        }
    }
}

#[allow(unused_macros)]
macro_rules! impl_neighbor_cartessian_nD{
    ($neighbor_name:ident, $cartessian_name:ident, $dim : expr) =>{
        doc_comment!{
            concat!(
                "Open-Ball neighborhood corresponds to [`", stringify!($cartessian_name), "`](struct@", stringify!($cartessian_name), ")"
            ),
            #[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
            pub struct $neighbor_name<'a>{
                /// Point of ball center
                pub center : &'a $cartessian_name<f64>,
            }
        }

        impl Neighbor<'_> for $neighbor_name<'_>{}

        impl<'a> $neighbor_name<'a>{
            doc_comment!{
                concat!(
                    "Initializing [`$neighbor_name`](struct@$neighbor_name) with a reference of the point.

# Examples

```
# use moldybrody::system::point::", stringify!($cartessian_name), ";
# use moldybrody::system::neighbor::", stringify!($neighbor_name), ";
let v = ", stringify!($cartessian_name), "{coord : [0.0; ", $dim, "]};
let n = ", stringify!($neighbor_name), "::new(&v);
```"
                ),
                pub fn new(center : &'a $cartessian_name<f64>) -> Self{
                    $neighbor_name{
                        center,
                    }
                }
            }
        }
    }
}


impl_neighbor_cartessian_nD!(OpenBall2D, Cartessian2D, 2);
impl_neighbor_cartessian_nD!(OpenBall3D, Cartessian3D, 3);
impl_neighbor_cartessian_nD!(OpenBall4D, Cartessian4D, 4);

// ==================================================================================
// ==================================================================================
// ==================================================================================


