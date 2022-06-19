use std::sync::Arc;

pub mod parallel;
pub mod serial;

pub trait PairIterator<'a, T: 'a>: Sized {
    type IntoIter: Iterator<Item = (&'a T, &'a T)>;
    type IntoEnum: Iterator<Item = ((usize, &'a T), (usize, &'a T))>;

    fn pair_iter(&'a self) -> Self::IntoIter;
    fn pair_enumerate(&'a self) -> Self::IntoEnum;
}

pub trait ParallelPairIterator<T>: Sized {
    type IntoIter: Iterator<Item = (Arc<T>, Arc<T>)>;
    type IntoEnum: Iterator<Item = ((usize, Arc<T>), (usize, Arc<T>))>;

    fn par_pair_iter(&self) -> Self::IntoIter;
    fn par_pair_enumerate(&self) -> Self::IntoEnum;
}
