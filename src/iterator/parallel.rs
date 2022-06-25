use std::collections::VecDeque;
use std::ops::Index;
use std::ops::Range;

use std::sync::{Arc, RwLock};

pub trait ParallelPairIterator<V>: Sized {
    type IntoIter: Iterator<Item = RowIter<V>>;
    type IntoEnum: Iterator<Item = RowEnum<V>>;

    fn par_pair_iter(&self) -> Self::IntoIter;
    fn par_pair_enumerate(&self) -> Self::IntoEnum;
}

impl<T> ParallelPairIterator<Vec<Arc<T>>> for Arc<Vec<Arc<T>>> {
    type IntoIter = IterDivider<Vec<Arc<T>>>;
    type IntoEnum = EnumDivider<Vec<Arc<T>>>;

    fn par_pair_iter(&self) -> IterDivider<Vec<Arc<T>>> {
        IterDivider::new(self, (*self).len())
    }

    fn par_pair_enumerate(&self) -> EnumDivider<Vec<Arc<T>>> {
        EnumDivider::new(self, (*self).len())
    }
}

impl<T: Copy, const N: usize> ParallelPairIterator<[Arc<T>; N]> for Arc<[Arc<T>; N]> {
    type IntoIter = IterDivider<[Arc<T>; N]>;
    type IntoEnum = EnumDivider<[Arc<T>; N]>;

    fn par_pair_iter(&self) -> IterDivider<[Arc<T>; N]> {
        IterDivider::new(self, self.len())
    }

    fn par_pair_enumerate(&self) -> EnumDivider<[Arc<T>; N]> {
        EnumDivider::new(self, self.len())
    }
}

impl<T> ParallelPairIterator<VecDeque<Arc<T>>> for Arc<VecDeque<Arc<T>>> {
    type IntoIter = IterDivider<VecDeque<Arc<T>>>;
    type IntoEnum = EnumDivider<VecDeque<Arc<T>>>;

    fn par_pair_iter(&self) -> IterDivider<VecDeque<Arc<T>>> {
        IterDivider::new(self, self.len())
    }

    fn par_pair_enumerate(&self) -> EnumDivider<VecDeque<Arc<T>>> {
        EnumDivider::new(self, self.len())
    }
}

#[derive(Clone, Debug)]
pub struct IterDivider<V> {
    items: Arc<V>,
    end: usize,
    row_iterator: Arc<RwLock<Range<usize>>>,
}

impl<V> IterDivider<V> {
    #[allow(dead_code)]
    fn new<'a>(items: &'a Arc<V>, end: usize) -> Self {
        Self {
            items: (*items).clone(),
            end,
            row_iterator: Arc::new(RwLock::new(0..end - 1)),
        }
    }

    pub fn full_len(&self) -> usize {
        let e = self.end;
        e * (e - 1) / 2
    }

    pub fn reset(&self) {
        let mut lock_iter = self.row_iterator.write().unwrap();
        lock_iter.start = 0;
    }

    pub fn kill(&self) {
        let mut lock_iter = self.row_iterator.write().unwrap();
        lock_iter.start = self.end - 1;
    }
}

impl<V> Iterator for IterDivider<V> {
    type Item = RowIter<V>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut lock_iter = self.row_iterator.write().unwrap();
        lock_iter
            .next()
            .map(|row| RowIter::new(&self.items, row, self.end))
    }
}

impl<V> ExactSizeIterator for IterDivider<V> {
    fn len(&self) -> usize {
        self.row_iterator.read().unwrap().len()
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct RowIter<V> {
    items: Arc<V>,
    row: usize,
    col_iterator: Range<usize>,
}

impl<V> RowIter<V> {
    fn new(items: &Arc<V>, row: usize, end: usize) -> Self {
        Self {
            items: (*items).clone(),
            row,
            col_iterator: (row + 1..end),
        }
    }
}

impl<V, T> Iterator for RowIter<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    type Item = (Arc<T>, Arc<T>);

    fn next(&mut self) -> Option<Self::Item> {
        self.col_iterator
            .next()
            .map(|i| (self.items[self.row].clone(), self.items[i].clone()))
    }
}

impl<V, T> ExactSizeIterator for RowIter<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    fn len(&self) -> usize {
        self.col_iterator.len()
    }
}

#[derive(Clone, Debug)]
pub struct EnumDivider<V> {
    items: Arc<V>,
    end: usize,
    row_iterator: Arc<RwLock<Range<usize>>>,
}

impl<V> EnumDivider<V> {
    #[allow(dead_code)]
    fn new<'a>(items: &'a Arc<V>, end: usize) -> Self {
        Self {
            items: (*items).clone(),
            end,
            row_iterator: Arc::new(RwLock::new(0..end - 1)),
        }
    }

    pub fn full_len(&self) -> usize {
        let e = self.end;
        e * (e - 1) / 2
    }

    pub fn reset(&self) {
        let mut lock_iter = self.row_iterator.write().unwrap();
        lock_iter.start = 0;
    }

    pub fn kill(&self) {
        let mut lock_iter = self.row_iterator.write().unwrap();
        lock_iter.start = self.end - 1;
    }
}

impl<V> Iterator for EnumDivider<V> {
    type Item = RowEnum<V>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut lock_iter = self.row_iterator.write().unwrap();
        lock_iter
            .next()
            .map(|row| RowEnum::new(&self.items, row, self.end))
    }
}

impl<V> ExactSizeIterator for EnumDivider<V> {
    fn len(&self) -> usize {
        self.row_iterator.read().unwrap().len()
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct RowEnum<V> {
    items: Arc<V>,
    row: usize,
    col_iterator: Range<usize>,
}

impl<V> RowEnum<V> {
    fn new(items: &Arc<V>, row: usize, end: usize) -> Self {
        Self {
            items: (*items).clone(),
            row,
            col_iterator: (row + 1..end),
        }
    }
}

impl<V, T> Iterator for RowEnum<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    type Item = ((usize, Arc<T>), (usize, Arc<T>));

    fn next(&mut self) -> Option<Self::Item> {
        self.col_iterator.next().map(|i| {
            (
                (self.row, self.items[self.row].clone()),
                (i, self.items[i].clone()),
            )
        })
    }
}

impl<V, T> ExactSizeIterator for RowEnum<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    fn len(&self) -> usize {
        self.col_iterator.len()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::sync::mpsc::channel;
    use std::thread;

    #[test]
    fn test_iter_serial() {
        let n = 4;
        let list: Arc<Vec<Arc<usize>>> = Arc::new((0usize..n).map(|x| Arc::new(x)).collect());
        let mut divider = list.par_pair_iter();

        let mut row_iter = divider.next().unwrap();
        assert_eq!(row_iter, RowIter::new(&list, 0, n));
        assert_eq!(row_iter.next(), Some((Arc::new(0), Arc::new(1))));
        assert_eq!(row_iter.next(), Some((Arc::new(0), Arc::new(2))));
        assert_eq!(row_iter.next(), Some((Arc::new(0), Arc::new(3))));
        assert_eq!(row_iter.next(), None);

        let mut row_iter = divider.next().unwrap();
        assert_eq!(row_iter, RowIter::new(&list, 1, n));
        assert_eq!(row_iter.next(), Some((Arc::new(1), Arc::new(2))));
        assert_eq!(row_iter.next(), Some((Arc::new(1), Arc::new(3))));
        assert_eq!(row_iter.next(), None);

        let mut row_iter = divider.next().unwrap();
        assert_eq!(row_iter, RowIter::new(&list, 2, n));
        assert_eq!(row_iter.next(), Some((Arc::new(2), Arc::new(3))));
        assert_eq!(row_iter.next(), None);

        assert_eq!(divider.next(), None);
    }

    #[test]
    fn test_enum_serial() {
        let n = 4;
        let list: Arc<Vec<Arc<usize>>> = Arc::new((0usize..n).map(|x| Arc::new(x)).collect());
        let mut divider = list.par_pair_enumerate();

        let mut row_enum = divider.next().unwrap();
        assert_eq!(row_enum, RowEnum::new(&list, 0, n));
        assert_eq!(row_enum.next(), Some(((0, Arc::new(0)), (1, Arc::new(1)))));
        assert_eq!(row_enum.next(), Some(((0, Arc::new(0)), (2, Arc::new(2)))));
        assert_eq!(row_enum.next(), Some(((0, Arc::new(0)), (3, Arc::new(3)))));
        assert_eq!(row_enum.next(), None);

        let mut row_enum = divider.next().unwrap();
        assert_eq!(row_enum, RowEnum::new(&list, 1, n));
        assert_eq!(row_enum.next(), Some(((1, Arc::new(1)), (2, Arc::new(2)))));
        assert_eq!(row_enum.next(), Some(((1, Arc::new(1)), (3, Arc::new(3)))));
        assert_eq!(row_enum.next(), None);

        let mut row_enum = divider.next().unwrap();
        assert_eq!(row_enum, RowEnum::new(&list, 2, n));
        assert_eq!(row_enum.next(), Some(((2, Arc::new(2)), (3, Arc::new(3)))));
        assert_eq!(row_enum.next(), None);

        assert_eq!(divider.next(), None);
    }

    #[test]
    fn test_parallel_iter_noflag() {
        let (n, n_thrd) = (5, 4);
        let (tx, rx) = channel();
        let mut threads = Vec::new();
        let mut results = Vec::new();

        let list: Arc<Vec<Arc<usize>>> = Arc::new((0usize..n).map(|x| Arc::new(x)).collect());
        let par_list = list.par_pair_iter();
        let mut remains = par_list.full_len();

        for _thread_num in 0..n_thrd {
            let thread_tx = tx.clone();
            let mut thread_list = par_list.clone();

            let handle = thread::spawn(move || {
                while let Some(iter) = thread_list.next() {
                    for (i, j) in iter {
                        thread_tx.send(*j - *i).unwrap();
                    }
                }
            });

            threads.push(handle);
        }

        while remains != 0 {
            match rx.recv() {
                Ok(x) => {
                    results.push(x);
                    remains -= 1;
                }
                Err(_) => {
                    break;
                }
            }
        }

        for handle in threads {
            handle.join().unwrap();
        }

        results.sort();
        assert_eq!(results, vec![1, 1, 1, 1, 2, 2, 2, 3, 3, 4]);
    }
}
