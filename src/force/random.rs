// macro_rules! impl_brownian {
//     ($ty : ident) => {
//         impl<'a, S, P> RandomForce<'a, $ty> for S
//         where
//             S: State<Position = P> + DiffusionFloat<$ty>,
//             P: Vector<Item = $ty> + Zeros + Map<Item = $ty>,
//             StandardNormal: Distribution<$ty>,
//         {
//             type Force = P;

//             fn force(&self, rng: &mut Pcg64) -> Self::Force {
//                 let mut force = P::zero_with_length(self.pos().dim());
//                 let c = ((2.0 as $ty) / self.diff_const()).sqrt();
//                 force.map_inplace(|x| *x = c * get_gaussian::<$ty>(rng));
//                 return force;
//             }

//             fn force_to(&self, rng: &mut Pcg64, force: &mut Self::Force) {
//                 let c = ((2.0 as $ty) / self.diff_const()).sqrt();
//                 force.map_inplace(|x| *x = c * get_gaussian::<$ty>(rng));
//             }

//             fn force_add_to(&self, rng: &mut Pcg64, force: &mut Self::Force) {
//                 let c = ((2.0 as $ty) / self.diff_const()).sqrt();
//                 force.map_inplace(|x| *x += c * get_gaussian::<$ty>(rng));
//             }
//         }

// impl<'a, S, P, const N : usize> RandomForce<'a, [$ty; N]> for S
//     where S : State<Position = P> + DiffusionFloat<[$ty; N]>,
//           P : Vector<Item = $ty> + Zeros + Map<Item = $ty>,
//           StandardNormal : Distribution<$ty>{
//     type Force = P;

//     fn force(&self, rng : &mut Pcg64) -> Self::Force {
//         let mut force = P::zero_with_length(N);
//         let diff = self.diff_const();
//         force.zip_mut_with(diff, |x, y| {
//             let c = ((2.0 as $ty) / y).sqrt();
//             *x = c * get_gaussian::<$ty>(rng);
//         });
//         return force;
//     }

//     fn force_to(&self, rng : &mut Pcg64, force : &mut Self::Force) {
//         let diff = self.diff_const();
//         force.zip_mut_with(diff, |x, y| {
//             let c = ((2.0 as $ty) / y).sqrt();
//             *x = c * get_gaussian::<$ty>(rng);
//         });
//     }

//     fn force_add_to(&self, rng : &mut Pcg64, force : &mut Self::Force) {
//         let diff = self.diff_const();
//         force.zip_mut_with(diff, |x, y| {
//             let c = ((2.0 as $ty) / y).sqrt();
//             *x += c * get_gaussian::<$ty>(rng);
//         });
//     }
// }
//     };
// }

// impl_brownian!(f32);
// impl_brownian!(f64);

// impl<'a, S, T, const N: usize> RandomForce<'a, ()> for S
// where
//     S: State<Position = Cartessian<T, N>>,
//     T: Scalar + Integer,
// {
//     type Force = [(f64, f64); N];

//     fn force(&self, _rng: &mut Pcg64) -> Self::Force {
//         let c: f64 = 1f64 / (2.0 * N as f64);
//         [(c, c); N]
//     }

//     fn force_to(&self, _rng: &mut Pcg64, force: &mut Self::Force) {
//         let c: f64 = 1f64 / (2.0 * N as f64);
//         force.iter_mut().for_each(|x| *x = (c, c));
//     }

//     fn force_add_to(&self, _rng: &mut Pcg64, force: &mut Self::Force) {
//         let c: f64 = 1f64 / (2.0 * N as f64);
//         force.iter_mut().for_each(|x| *x = (c, c));
//     }
// }
