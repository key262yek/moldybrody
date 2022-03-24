//! Module for random number generator
//!
//! 이 crate에는 여러 random number generator가 쓰입니다.
//! - Initial configuration : Uniform distribution
//! - Random Force : White noise, Levy walk
//! - Boundary condition : Partially absorbing boundary, Scattering boundary
//!
//! 각각의 case에서 모두 같은 [rand crate](https://docs.rs/rand/0.8.4/rand/)의 Pcg64를 random number generator로 씁니다.
//! 이를 declare하고 각 module에서 사용할 수 있도록 wrapping 해줍니다.

use rand::prelude::Distribution;
use crate::prelude::*;
use rand_distr::{StandardNormal, Open01};
use rand::Rng;



/// Declare Pcg64 random number generator from a given seed
///
/// PCG family algorithm을 기반해서 random number를 만드는 generator를 만들어주는 함수.
///
/// # Examples
/// ```
/// # use moldybrody::prelude::*;
/// use moldybrody::rng::rng_seed;
/// let seed : u128 = 12412398120480;
/// let rng = rng_seed(seed);
/// ```
pub fn rng_seed(seed : u128) -> Pcg64{
    const INC: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;
    rand_pcg::Pcg64::new(seed, INC)
}

/// Generate uniform distriubted random number
///
/// (0, 1) 범위의 random number를 출력해주는 함수
///
/// # Examples
/// ```
/// # use moldybrody::prelude::*;
/// # use moldybrody::rng::rng_seed;
/// use moldybrody::rng::get_uniform;
///
/// let seed : u128 = 12412398120480;
/// let mut rng = rng_seed(seed);
/// let uni = get_uniform(&mut rng);
/// ```
pub fn get_uniform<T>(rng : &mut Pcg64) -> T
     where Open01: Distribution<T>{
    rng.sample(Open01)
}

/// Generate gaussian random number
///
/// standard normal random number를 만들어주는 함수
///
/// # Examples
/// ```
/// # use moldybrody::prelude::*;
/// # use moldybrody::rng::rng_seed;
/// use moldybrody::rng::get_gaussian;
///
/// let seed : u128 = 12412398120480;
/// let mut rng = rng_seed(seed);
/// let gas = get_gaussian(&mut rng);
/// ```
pub fn get_gaussian<T>(rng : &mut Pcg64) -> T
    where StandardNormal : Distribution<T>{
    rng.sample(StandardNormal)
}
