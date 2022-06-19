// use rayon::prelude::IntoParallelIterator;
use super::PairIterator;
use std::collections::VecDeque;

use std::ops::Index;
use std::ops::Range;

#[derive(Clone)]
pub struct PairIdx {
    start: usize,
    end: usize,
    current: (usize, usize),
}

impl PairIdx {
    #[allow(dead_code)]
    fn new(start: usize, end: usize) -> Self {
        Self {
            start,
            end,
            current: (start, start + 1),
        }
    }

    pub fn reset(&mut self) {
        let s = self.start;
        self.current = (s, s + 1);
    }

    pub fn kill(&mut self) {
        let e = self.end;
        self.current = (e - 1, e);
    }

    pub fn full_len(&self) -> usize {
        let e = self.end - self.start;
        e * (e - 1) / 2
    }
}

impl Iterator for PairIdx {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current.1 == self.end {
            return None;
        } else if self.current.1 == self.end - 1 {
            let i = self.current.0;
            self.current = (i + 1, i + 2);
            return Some((i, self.end - 1));
        } else {
            let tmp = self.current;
            self.current.1 += 1;
            return Some(tmp);
        }
    }
}

impl ExactSizeIterator for PairIdx {
    fn len(&self) -> usize {
        let (e, c0, c1) = (
            self.end as i32,
            self.current.0 as i32,
            self.current.1 as i32,
        );
        ((e - c0) * (e - c0 - 1) / 2 - c1 + c0 + 1) as usize
    }

    // fn is_empty(&self) -> bool {
    //     if self.current == (self.end, self.end + 1) {
    //         return true;
    //     } else {
    //         return false;
    //     }
    // }
}

pub trait PairIndex {
    type IntoIter: Iterator<Item = (usize, usize)>;

    fn pair(&self) -> Self::IntoIter;
}

impl<T> PairIndex for Range<T>
where
    T: Into<usize> + Copy,
{
    type IntoIter = PairIdx;

    fn pair(&self) -> PairIdx {
        PairIdx::new(self.start.into(), self.end.into())
    }
}

// impl<'a, T, V> PairIterator<'a, T> for V
// where
//     V: ExactSizeIterator + Index<usize, Output = T> + Sized,
// {
//     fn pair_iter(&'a self) -> PairIter<'a, Self> {
//         PairIter::new(self, self.len())
//     }

//     fn pair_enumerate(&'a self) -> PairEnum<'a, Self> {
//         PairEnum::new(self, self.len())
//     }
// }

#[derive(Clone)]
pub struct PairIter<'a, V> {
    items: &'a V,
    idx_iterator: PairIdx,
}

impl<'a, V> PairIter<'a, V> {
    fn new(items: &'a V, length: usize) -> Self {
        Self {
            items,
            idx_iterator: (0..length).pair(),
        }
    }

    pub fn reset(&mut self) {
        self.idx_iterator.reset();
    }

    pub fn kill(&mut self) {
        self.idx_iterator.kill();
    }

    pub fn full_len(&mut self) -> usize {
        self.idx_iterator.full_len()
    }
}

impl<'a, V> Iterator for PairIter<'a, V>
where
    V: Index<usize>,
    <V as Index<usize>>::Output: Sized,
{
    type Item = (
        &'a <V as Index<usize>>::Output,
        &'a <V as Index<usize>>::Output,
    );

    fn next(&mut self) -> Option<Self::Item> {
        match self.idx_iterator.next() {
            Some((i, j)) => {
                return Some((&self.items[i], &self.items[j]));
            }
            None => {
                return None;
            }
        }
    }
}

impl<'a, V> ExactSizeIterator for PairIter<'a, V>
where
    V: Index<usize>,
    <V as Index<usize>>::Output: Sized,
{
    fn len(&self) -> usize {
        self.idx_iterator.len()
    }
}

// impl<'a, V> IntoParallelIterator for PairIter<'a, V>
// where
//     V: Index<usize>,
//     <V as Index<usize>>::Output: Sized + Send,
// {
//     type Item = (
//         &'a <V as Index<usize>>::Output,
//         &'a <V as Index<usize>>::Output,
//     );
//     type Iter = rayon::vec::Iter<'a, ;

//     fn into_par_iter(self) -> Self::Iter {
//         self.iter()
//             .collect::<Vec<(
//                 &'a <V as Index<usize>>::Output,
//                 &'a <V as Index<usize>>::Output,
//             )>>()
//             .into_par_iter()
//     }
// }

#[derive(Clone)]
pub struct PairEnum<'a, V> {
    items: &'a V,
    idx_iterator: PairIdx,
}

impl<'a, V> PairEnum<'a, V> {
    fn new(items: &'a V, length: usize) -> Self {
        Self {
            items,
            idx_iterator: (0..length).pair(),
        }
    }

    pub fn reset(&mut self) {
        self.idx_iterator.reset();
    }

    pub fn kill(&mut self) {
        self.idx_iterator.kill();
    }

    pub fn full_len(&mut self) -> usize {
        self.idx_iterator.full_len()
    }
}

impl<'a, V> Iterator for PairEnum<'a, V>
where
    V: Index<usize>,
    <V as Index<usize>>::Output: Sized,
{
    type Item = (
        (usize, &'a <V as Index<usize>>::Output),
        (usize, &'a <V as Index<usize>>::Output),
    );

    fn next(&mut self) -> Option<Self::Item> {
        match self.idx_iterator.next() {
            Some((i, j)) => {
                return Some(((i, &self.items[i]), (j, &self.items[j])));
            }
            None => {
                return None;
            }
        }
    }
}

impl<'a, V> ExactSizeIterator for PairEnum<'a, V>
where
    V: Index<usize>,
    <V as Index<usize>>::Output: Sized,
{
    fn len(&self) -> usize {
        self.idx_iterator.len()
    }
}

impl<'a, T: 'a> PairIterator<'a, T> for Vec<T> {
    type IntoIter = PairIter<'a, Self>;
    type IntoEnum = PairEnum<'a, Self>;

    fn pair_iter(&'a self) -> PairIter<'a, Self> {
        PairIter::new(self, self.len())
    }

    fn pair_enumerate(&'a self) -> PairEnum<'a, Self> {
        PairEnum::new(self, self.len())
    }
}

impl<'a, T: 'a, const N: usize> PairIterator<'a, T> for [T; N] {
    type IntoIter = PairIter<'a, Self>;
    type IntoEnum = PairEnum<'a, Self>;

    fn pair_iter(&'a self) -> PairIter<'a, Self> {
        PairIter::new(self, self.len())
    }

    fn pair_enumerate(&'a self) -> PairEnum<'a, Self> {
        PairEnum::new(self, self.len())
    }
}

impl<'a, T: 'a> PairIterator<'a, T> for VecDeque<T> {
    type IntoIter = PairIter<'a, Self>;
    type IntoEnum = PairEnum<'a, Self>;

    fn pair_iter(&'a self) -> PairIter<'a, Self> {
        PairIter::new(self, self.len())
    }

    fn pair_enumerate(&'a self) -> PairEnum<'a, Self> {
        PairEnum::new(self, self.len())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rayon::prelude::*;

    #[test]
    fn test_pair() {
        let mut iter: PairIdx = (3usize..6).pair();
        assert_eq!(iter.next(), Some((3, 4)));
        assert_eq!(iter.next(), Some((3, 5)));
        assert_eq!(iter.next(), Some((4, 5)));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_vec_pair_iter() {
        let array = [3, 4, 5, 6];
        let mut iter = array.pair_iter();
        assert_eq!(iter.next(), Some((&3, &4)));
        assert_eq!(iter.next(), Some((&3, &5)));
        assert_eq!(iter.next(), Some((&3, &6)));
        assert_eq!(iter.next(), Some((&4, &5)));
        assert_eq!(iter.next(), Some((&4, &6)));
        assert_eq!(iter.next(), Some((&5, &6)));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_vec_pair_enum() {
        let array = [3, 4, 5, 6];
        let mut iter = array.pair_enumerate();
        assert_eq!(iter.next(), Some(((0, &3), (1, &4))));
        assert_eq!(iter.next(), Some(((0, &3), (2, &5))));
        assert_eq!(iter.next(), Some(((0, &3), (3, &6))));
        assert_eq!(iter.next(), Some(((1, &4), (2, &5))));
        assert_eq!(iter.next(), Some(((1, &4), (3, &6))));
        assert_eq!(iter.next(), Some(((2, &5), (3, &6))));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_pair_par() {
        let array = [3, 4, 5, 6];
        let mut items: Vec<(&i32, &i32)> = array.pair_iter().par_bridge().collect();
        items.par_sort();
        assert_eq!(
            items,
            vec![(&3, &4), (&3, &5), (&3, &6), (&4, &5), (&4, &6), (&5, &6)]
        );
    }

    #[test]
    fn test_exact_size() {
        let mut iter: PairIdx = (3usize..6).pair();
        assert_eq!(iter.len(), 3);
        iter.next();
        assert_eq!(iter.len(), 2);
        iter.next();
        assert_eq!(iter.len(), 1);
        iter.next();
        assert_eq!(iter.len(), 0);
        iter.next();
        assert_eq!(iter.len(), 0);
    }

    #[test]
    fn test_reset_kill() {
        let mut iter: PairIdx = (3usize..6).pair();
        assert_eq!(iter.next(), Some((3, 4)));
        iter.kill();
        assert_eq!(iter.next(), None);

        iter.reset();
        assert_eq!(iter.next(), Some((3, 4)));
    }
}
