use std::collections::VecDeque;
use std::ops::Index;
use std::ops::Range;

use std::sync::{Arc, RwLock};

trait ParallelPairIterator<T>: Sized {
    type IntoIter: Iterator<Item = (Arc<T>, Arc<T>)>;
    type IntoEnum: Iterator<Item = ((usize, Arc<T>), (usize, Arc<T>))>;

    fn par_pair_iter(&self) -> Self::IntoIter;
    fn par_pair_enumerate(&self) -> Self::IntoEnum;
}

pub trait ParallelPairIndex {
    type IntoIter: Iterator<Item = (usize, usize)>;

    fn par_pair(&self) -> Self::IntoIter;
}

#[derive(Clone)]
pub struct ParallelPairIdx {
    start: usize,
    end: usize,
    current: Arc<RwLock<(usize, usize)>>,
}

impl ParallelPairIdx {
    #[allow(dead_code)]
    fn new(start: usize, end: usize) -> Self {
        Self {
            start,
            end,
            current: Arc::new(RwLock::new((start, start + 1))),
        }
    }

    pub fn reset(&self) {
        let s = self.start;
        let mut lock_current = self.current.write().unwrap();
        *lock_current = (s, s + 1);
    }

    pub fn kill(&self) {
        let e = self.end;
        let mut lock_current = self.current.write().unwrap();
        *lock_current = (e - 1, e);
    }

    pub fn full_len(&self) -> usize {
        let e = self.end - self.start;
        e * (e - 1) / 2
    }
}

impl ExactSizeIterator for ParallelPairIdx {
    fn len(&self) -> usize {
        let c = self.current.read().unwrap();
        let (e, c0, c1) = (self.end as i32, c.0 as i32, c.1 as i32);
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

impl Iterator for ParallelPairIdx {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let mut mutex_current = self.current.write().unwrap();

        if mutex_current.1 == self.end {
            return None;
        } else if mutex_current.1 == self.end - 1 {
            let i = mutex_current.0;
            *mutex_current = (i + 1, i + 2);
            return Some((i, self.end - 1));
        } else {
            let tmp = *mutex_current;
            mutex_current.1 += 1;
            return Some(tmp);
        }
    }
}

impl<T> ParallelPairIndex for Range<T>
where
    T: Into<usize> + Copy,
{
    type IntoIter = ParallelPairIdx;

    fn par_pair(&self) -> ParallelPairIdx {
        ParallelPairIdx::new(self.start.into(), self.end.into())
    }
}

#[derive(Clone)]
pub struct ParallelPairIter<V> {
    items: Arc<V>,
    idx_iterator: Arc<RwLock<ParallelPairIdx>>,
}

impl<V> ParallelPairIter<V> {
    fn new<'a>(items: &'a Arc<V>, length: usize) -> Self {
        Self {
            items: (*items).clone(),
            idx_iterator: Arc::new(RwLock::new((0..length).par_pair())),
        }
    }

    pub fn reset(&self) {
        self.idx_iterator.write().unwrap().reset();
    }

    pub fn kill(&self) {
        self.idx_iterator.write().unwrap().kill();
    }

    pub fn full_len(&self) -> usize {
        self.idx_iterator.read().unwrap().full_len()
    }
}

impl<V, T> Iterator for ParallelPairIter<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    type Item = (Arc<T>, Arc<T>);

    fn next(&mut self) -> Option<Self::Item> {
        let mut iter = self.idx_iterator.write().unwrap();
        match iter.next() {
            Some((i, j)) => {
                return Some((self.items[i].clone(), self.items[j].clone()));
            }
            None => {
                return None;
            }
        }
    }
}

impl<V, T> ExactSizeIterator for ParallelPairIter<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    fn len(&self) -> usize {
        let iter = self.idx_iterator.read().unwrap();
        iter.len()
    }
}

#[derive(Clone)]
pub struct ParallelPairEnum<V> {
    items: Arc<V>,
    idx_iterator: Arc<RwLock<ParallelPairIdx>>,
}

impl<V> ParallelPairEnum<V> {
    fn new<'a>(items: &'a Arc<V>, length: usize) -> Self {
        Self {
            items: (*items).clone(),
            idx_iterator: Arc::new(RwLock::new((0..length).par_pair())),
        }
    }

    pub fn reset(&self) {
        self.idx_iterator.write().unwrap().reset();
    }

    pub fn kill(&self) {
        self.idx_iterator.write().unwrap().kill();
    }

    pub fn full_len(&self) -> usize {
        self.idx_iterator.read().unwrap().full_len()
    }
}

impl<V, T> Iterator for ParallelPairEnum<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    type Item = ((usize, Arc<T>), (usize, Arc<T>));

    fn next(&mut self) -> Option<Self::Item> {
        let mut iter = self.idx_iterator.write().unwrap();
        match iter.next() {
            Some((i, j)) => {
                return Some(((i, self.items[i].clone()), (j, self.items[j].clone())));
            }
            None => {
                return None;
            }
        }
    }
}

impl<V, T> ExactSizeIterator for ParallelPairEnum<V>
where
    V: Index<usize, Output = Arc<T>>,
{
    fn len(&self) -> usize {
        let iter = self.idx_iterator.read().unwrap();
        iter.len()
    }
}

// impl<T> ParallelPairIterator<T> for Arc<Vec<Arc<T>>> {
//     type IntoIter = ParallelPairIter<Vec<Arc<T>>>;
//     type IntoEnum = ParallelPairEnum<Vec<Arc<T>>>;

//     fn par_pair_iter(&self) -> ParallelPairIter<Vec<Arc<T>>> {
//         ParallelPairIter::new(self, (*self).len())
//     }

//     fn par_pair_enumerate(&self) -> ParallelPairEnum<Vec<Arc<T>>> {
//         ParallelPairEnum::new(self, (*self).len())
//     }
// }

impl<T: Copy, const N: usize> ParallelPairIterator<T> for Arc<[Arc<T>; N]> {
    type IntoIter = ParallelPairIter<[Arc<T>; N]>;
    type IntoEnum = ParallelPairEnum<[Arc<T>; N]>;

    fn par_pair_iter(&self) -> ParallelPairIter<[Arc<T>; N]> {
        ParallelPairIter::new(self, self.len())
    }

    fn par_pair_enumerate(&self) -> ParallelPairEnum<[Arc<T>; N]> {
        ParallelPairEnum::new(self, self.len())
    }
}

impl<T> ParallelPairIterator<T> for Arc<VecDeque<Arc<T>>> {
    type IntoIter = ParallelPairIter<VecDeque<Arc<T>>>;
    type IntoEnum = ParallelPairEnum<VecDeque<Arc<T>>>;

    fn par_pair_iter(&self) -> ParallelPairIter<VecDeque<Arc<T>>> {
        ParallelPairIter::new(self, self.len())
    }

    fn par_pair_enumerate(&self) -> ParallelPairEnum<VecDeque<Arc<T>>> {
        ParallelPairEnum::new(self, self.len())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::sync::mpsc::channel;
    use std::sync::Mutex;
    use std::thread;

    #[test]
    fn test_parallel_index() {
        fn check(flag: &Arc<Mutex<bool>>) -> bool {
            let b = flag.lock().unwrap();
            *b
        }
        let (n, n_thrd) = (10, 4);
        let (tx, rx) = channel();
        let mut threads = Vec::new();
        let mut results = Vec::new();
        let mut remains = n * (n - 1) / 2;
        let flag = Arc::new(Mutex::new(true));

        let par_pair = (0..n).par_pair();

        for _thread_num in 0..n_thrd {
            let thread_tx = tx.clone();
            let mut thread_queue = par_pair.clone();
            let thread_flag = flag.clone();

            let handle = thread::spawn(move || {
                while check(&thread_flag) {
                    // println!("{:?}", thread_queue.current);
                    if let Some((i, j)) = thread_queue.next() {
                        thread_tx.send((i, j)).unwrap();
                    } else {
                        break;
                    }
                }
                // println!("While loop end!");
            });

            threads.push(handle);
        }

        while remains != 0 {
            match rx.recv() {
                Ok(x) => {
                    results.push(x);
                    remains -= 1;
                    // println!("{:?}", remains);
                }
                Err(_) => {
                    break;
                }
            }
        }

        {
            let mut b = flag.lock().unwrap();
            *b = false;
        }

        for handle in threads {
            handle.join().unwrap();
        }

        // println!("{:?}", results);
        results.sort();
        assert_eq!(
            results,
            (0usize..n).par_pair().collect::<Vec<(usize, usize)>>()
        );
    }

    #[test]
    fn test_parallel_index_noflag() {
        let (n, n_thrd) = (10, 4);
        let (tx, rx) = channel();
        let mut threads = Vec::new();
        let mut results = Vec::new();
        let mut remains = n * (n - 1) / 2;

        let par_pair = (0..n).par_pair();

        for _thread_num in 0..n_thrd {
            let thread_tx = tx.clone();
            let mut thread_queue = par_pair.clone();

            let handle = thread::spawn(move || {
                while let Some((i, j)) = thread_queue.next() {
                    thread_tx.send((i, j)).unwrap();
                }
                // println!("While loop end!");
            });

            threads.push(handle);
        }

        while remains != 0 {
            match rx.recv() {
                Ok(x) => {
                    results.push(x);
                    remains -= 1;
                    // println!("{:?}", remains);
                }
                Err(_) => {
                    break;
                }
            }
        }

        for handle in threads {
            handle.join().unwrap();
        }

        // println!("{:?}", results);
        results.sort();
        assert_eq!(
            results,
            (0usize..n).par_pair().collect::<Vec<(usize, usize)>>()
        );
    }

    #[test]
    fn test_parallel_exact_size() {
        let mut iter: ParallelPairIdx = (3usize..6).par_pair();
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
    fn test_parallel_iter_noflag() {
        let (n, n_thrd) = (5, 4);
        let (tx, rx) = channel();
        let mut threads = Vec::new();
        let mut results = Vec::new();

        let list: Arc<Vec<Arc<usize>>> = Arc::new((0usize..n).map(|x| Arc::new(x)).collect());
        let par_list = ParallelPairIter::new(&list, n);
        let mut remains = par_list.len();

        for _thread_num in 0..n_thrd {
            let thread_tx = tx.clone();
            let mut thread_list = par_list.clone();

            let handle = thread::spawn(move || {
                while let Some((i, j)) = thread_list.next() {
                    thread_tx.send(*j - *i).unwrap();
                }
                // println!("While loop end!");
            });

            threads.push(handle);
        }

        while remains != 0 {
            match rx.recv() {
                Ok(x) => {
                    results.push(x);
                    remains -= 1;
                    // println!("{:?}", remains);
                }
                Err(_) => {
                    break;
                }
            }
        }

        for handle in threads {
            handle.join().unwrap();
        }

        // println!("{:?}", results);
        results.sort();
        assert_eq!(results, vec![1, 1, 1, 1, 2, 2, 2, 3, 3, 4]);
    }
}
