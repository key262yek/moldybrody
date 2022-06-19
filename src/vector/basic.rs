

// use std::convert::TryInto;
// use std::fmt::Debug;
use num_traits::Zero;
use crate::vector::{Cartessian, CartessianND};

use std::{
    slice::{Iter, IterMut},
    ops::{Index, IndexMut, FnMut},
    hash::{Hash, Hasher},
    borrow::{Borrow, BorrowMut},
};

pub trait Zeros{
    fn zero_with_length(n : usize) -> Self;
}

pub trait Map{
    type Item : Clone;

    fn map_inplace<'a, F>(&'a mut self, f : F)
        where F : FnMut(&'a mut Self::Item);

    fn zip_mut_with<F, U>(&mut self, rhs : U, f : F)
    where F : FnMut(&mut Self::Item, <U as IntoIterator>::Item),
          U : IntoIterator;
}

impl<T, const N : usize> Cartessian<T, N>{
    /// Return a reference to the element of coordinate at index, or return None if the index is out of bounds.
    /// Arrays also support indexing syntax: array[index].
    ///
    /// ```
    /// # use moldybrody::vector::Cartessian;
    /// let a = Cartessian::<f64, 4>::new([2.0, 3.0, 4.0, 5.0]);
    ///
    /// assert!(
    ///     a.get(0) == Some(&2.0) &&
    ///     a.get(3) == Some(&5.0) &&
    ///     a[1] == 3.0 &&
    ///     a[2] == 4.0
    /// )
    /// ```
    pub fn get(&self, index : usize) -> Option<&T>{
        if index >= N {
            return None;
        }
        return Some(&self.coord[index])
    }

    /// Return a mutable reference to the element at index, or return None if the index is out of bounds.
    pub fn get_mut(&mut self, index : usize) -> Option<&mut T>{
        if index >= N {
            return None;
        }
        return Some(&mut self.coord[index])
    }

    // Call f by reference on each element and create a new vector with the new values.
    // Return an array with the same shape as self.
    //
    // ```
    // # use moldybrody::vector::Cartessian;
    // let a : Cartessian<f64, 2> = Cartessian::new([2.0, 0.0]);
    // assert!(
    //     a.map(|x| *x as i32) == Cartessian::<i32, 2>::new([2, 0])
    // )
    // ```
    // pub fn map<'a, B, F>(&'a self, mut f : F) -> Cartessian<B, N>
    //     where F : FnMut(&'a T) -> B,
    //           T : Clone,
    //           B : Debug{

    //     // Cartessian::<B, N>{
    //     //     coord : self.into_iter().map(|x| f(x)).collect::<Vec<B>>().try_into().unwrap(),
    //     // }

    //     Cartessian::<B, N>::from_iter_unchecked(self.into_iter().map(|x| f(x)))
    // }

    // Modify the array in place by calling f by mutable reference on each element.
    // pub fn map_inplace<'a, F>(&'a mut self, f : F)
    //     where F : FnMut(&'a mut T){

    //     self.iter_mut().for_each(f)
    // }

    // /// calling the closure f on each element pair.
    // pub fn zip_mut_with<F>(&mut self, rhs : &Self, mut f : F)
    //     where T : Clone,
    //           F : FnMut(&mut T, &T){

    //     self.iter_mut().zip(rhs).for_each(|(x, y)| f(x, y));
    // }
}

impl<T, const N : usize> Map for Cartessian<T, N>
    where T : Clone{
    type Item = T;

    fn map_inplace<'a, F>(&'a mut self, f : F)
        where F : FnMut(&'a mut T) {
        self.iter_mut().for_each(f)
    }

    fn zip_mut_with<F, U>(&mut self, rhs : U, mut f : F)
        where F : FnMut(&mut T, <U as IntoIterator>::Item),
            U : IntoIterator {
        self.iter_mut().zip(rhs).for_each(|(x, y)| f(x, y));
    }
}


impl<T, const N : usize> AsMut<[T]> for Cartessian<T, N>{
    fn as_mut(&mut self) -> &mut [T]{
        self.coord.as_mut()
    }
}

impl<T, const N : usize> AsRef<[T]> for Cartessian<T, N>{
    fn as_ref(&self) -> &[T]{
        self.coord.as_ref()
    }
}

impl<T, const N : usize> Borrow<[T]> for Cartessian<T, N>{
    fn borrow(&self) -> &[T]{
        self.coord.borrow()
    }
}

impl<T, const N : usize> BorrowMut<[T]> for Cartessian<T, N>{
    fn borrow_mut(&mut self) -> &mut [T]{
        self.coord.borrow_mut()
    }
}

impl<T, const N : usize> Default for Cartessian<T, N>
    where T : Default + Copy{
    fn default() -> Self{
        Self{
            coord : [T::default(); N],
        }
    }
}

impl<T, const N : usize> Zeros for Cartessian<T, N>
    where T : Zero + Copy{
    fn zero_with_length(n : usize) -> Self {
        if n != N{
            panic!("Invalid Argument : Cartessian<T, {}> cannot be defined with length {}", N, n);
        }
        Self{
            coord : [T::zero(); N],
        }
    }
}

impl<T, const N : usize> Hash for Cartessian<T, N>
    where T : Hash{
    fn hash<H>(&self, state : &mut H) where H : Hasher{
        self.coord.hash(state)
    }
}

impl<T, I, const N : usize> Index<I> for Cartessian<T, N>
    where [T]:Index<I>{
    type Output = <[T] as Index<I>>::Output;

    fn index(&self, index : I) -> &<[T] as Index<I>>::Output{
        self.coord.index(index)
    }
}

impl<T, I, const N : usize> IndexMut<I> for Cartessian<T, N>
    where [T]:IndexMut<I>{

    fn index_mut(&mut self, index : I) -> &mut <[T] as Index<I>>::Output{
        self.coord.index_mut(index)
    }
}

impl<T, const N : usize> Cartessian<T, N>{
    pub fn iter<'a>(&'a self) -> Iter<'a, T>{
        self.coord.iter()
    }

    pub fn iter_mut<'a>(&'a mut self) -> IterMut<'a, T>{
        self.coord.iter_mut()
    }
}


impl<'a, T, const N : usize> IntoIterator for &'a Cartessian<T, N>{
    type Item = &'a T;
    type IntoIter = Iter<'a, T>;

    fn into_iter(self) -> Iter<'a, T>{
        self.coord.iter()
    }
}

impl<'a, T, const N : usize> IntoIterator for &'a mut Cartessian<T, N>{
    type Item = &'a mut T;
    type IntoIter = IterMut<'a, T>;

    fn into_iter(self) -> IterMut<'a, T>{
        self.coord.iter_mut()
    }
}


// ===========================================================================================================
// ===========================================================================================================

impl<T> CartessianND<T>{
    /// Return a reference to the element of coordinate at index, or return None if the index is out of bounds.
    /// Arrays also support indexing syntax: array[index].
    ///
    /// ```
    /// # use moldybrody::vector::CartessianND;
    /// let a = CartessianND::<f64>::new(vec![2.0, 3.0, 4.0, 5.0]);
    ///
    /// assert!(
    ///     a.get(0) == Some(&2.0) &&
    ///     a.get(3) == Some(&5.0) &&
    ///     a[1] == 3.0 &&
    ///     a[2] == 4.0
    /// )
    /// ```
    pub fn get(&self, index : usize) -> Option<&T>{
        if index >= self.len() {
            return None;
        }
        return Some(&self.coord[index])
    }

    /// Return a mutable reference to the element at index, or return None if the index is out of bounds.
    pub fn get_mut(&mut self, index : usize) -> Option<&mut T>{
        if index >= self.len() {
            return None;
        }
        return Some(&mut self.coord[index])
    }

    // Call f by reference on each element and create a new vector with the new values.
    // Return an array with the same shape as self.
    //
    // ```
    // # use moldybrody::vector::CartessianND;
    // let a : CartessianND<f64> = CartessianND::new(vec![2.0, 0.0]);
    // assert!(
    //     a.map(|x| *x as i32) == CartessianND::<i32>::new(vec![2, 0])
    // )
    // ```
    // pub fn map<'a, B, F>(&'a self, mut f : F) -> CartessianND<B>
    //     where F : FnMut(&'a T) -> B,
    //           B : Debug + Clone{

    //     CartessianND::<B>::from_iter(self.into_iter().map(|x| f(x)))
    // }

    // Modify the array in place by calling f by mutable reference on each element.
    // pub fn map_inplace<'a, F>(&'a mut self, f : F)
    //     where F : FnMut(&'a mut T){

    //     self.iter_mut().for_each(f)
    // }

    // /// calling the closure f on each element pair.
    // pub fn zip_mut_with<F>(&mut self, rhs : &Self, mut f : F)
    //     where T : Clone,
    //           F : FnMut(&mut T, &T){

    //     self.iter_mut().zip(rhs).for_each(|(x, y)| f(x, y));
    // }
}

impl<T> Map for CartessianND<T>
    where T : Clone{
    type Item = T;

    fn map_inplace<'a, F>(&'a mut self, f : F)
        where F : FnMut(&'a mut T) {
        self.iter_mut().for_each(f)
    }

    fn zip_mut_with<F, U>(&mut self, rhs : U, mut f : F)
        where F : FnMut(&mut T, <U as IntoIterator>::Item),
        U : IntoIterator {
        self.iter_mut().zip(rhs).for_each(|(x, y)| f(x, y));
    }
}


impl<T> AsMut<Vec<T>> for CartessianND<T>{
    fn as_mut(&mut self) -> &mut Vec<T>{
        self.coord.as_mut()
    }
}

impl<T> AsRef<Vec<T>> for CartessianND<T>{
    fn as_ref(&self) -> &Vec<T>{
        self.coord.as_ref()
    }
}

impl<T> Borrow<Vec<T>> for CartessianND<T>{
    fn borrow(&self) -> &Vec<T>{
        self.coord.borrow()
    }
}

impl<T> BorrowMut<Vec<T>> for CartessianND<T>{
    fn borrow_mut(&mut self) -> &mut Vec<T>{
        self.coord.borrow_mut()
    }
}


impl<T> Hash for CartessianND<T>
    where T : Hash{
    fn hash<H>(&self, state : &mut H) where H : Hasher{
        self.coord.hash(state)
    }
}

impl<T, I> Index<I> for CartessianND<T>
    where Vec<T> : Index<I>{
    type Output = <Vec<T> as Index<I>>::Output;

    fn index(&self, index : I) -> &<Vec<T> as Index<I>>::Output{
        self.coord.index(index)
    }
}

impl<T, I> IndexMut<I> for CartessianND<T>
    where Vec<T> : IndexMut<I>{

    fn index_mut(&mut self, index : I) -> &mut <Vec<T>  as Index<I>>::Output{
        self.coord.index_mut(index)
    }
}

impl<T> CartessianND<T>{
    pub fn iter<'a>(&'a self) -> Iter<'a, T>{
        self.coord.iter()
    }

    pub fn iter_mut<'a>(&'a mut self) -> IterMut<'a, T>{
        self.coord.iter_mut()
    }
}


impl<'a, T> IntoIterator for &'a CartessianND<T>{
    type Item = &'a T;
    type IntoIter = Iter<'a, T>;

    fn into_iter(self) -> Iter<'a, T>{
        self.coord.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut CartessianND<T>{
    type Item = &'a mut T;
    type IntoIter = IterMut<'a, T>;

    fn into_iter(self) -> IterMut<'a, T>{
        self.coord.iter_mut()
    }
}

impl<T> Zeros for CartessianND<T>
    where T : Zero + Copy{
    fn zero_with_length(n : usize) -> Self {
        CartessianND::<T>::new(vec![T::zero(); n])
    }
}

#[cfg(test)]
mod test {
    use crate::vector::basic::Map;
use crate::vector::{CartessianND, Cartessian2D};

    #[test]
    fn test_basic(){
        let mut a = Cartessian2D::new([1.0f64, 3.0]);
        assert_eq!(a.get(0), Some(&1.0));
        assert_eq!(a.get(1), Some(&3.0));

        *a.get_mut(0).unwrap() = 2.0;
        assert_eq!(a.get(0), Some(&2.0));

        // assert_eq!(a.map(|x| x * x), Cartessian2D::new([4.0, 9.0]));

        a.map_inplace(|x| *x = *x * *x);
        assert_eq!(a, Cartessian2D::new([4.0, 9.0]));

        let b = Cartessian2D::new([1.0, 3.0]);
        a.zip_mut_with(&b, |x, y| {*x = *x + *y});
        assert_eq!(a, Cartessian2D::new([5.0, 12.0]));

        assert_eq!(a.as_ref(), &[5.0, 12.0]);

        let a_mut = a.as_mut();
        a_mut[0] = 1.0;
        assert_eq!(a, Cartessian2D::new([1.0, 12.0]));

        assert_eq!(a[0], 1.0);
        assert_eq!(a[1], 12.0);
        a[0] = 3.0;
        a[1] = 9.0;
        assert_eq!(a, Cartessian2D::new([3.0, 9.0]));

        let mut iter = a.into_iter();
        assert_eq!(iter.next(), Some(&3.0));
        assert_eq!(iter.next(), Some(&9.0));
        assert_eq!(iter.next(), None);




        let mut c = CartessianND::new(vec![1.0f64, 3.0]);
        assert_eq!(c.get(0), Some(&1.0));
        assert_eq!(c.get(1), Some(&3.0));

        *c.get_mut(0).unwrap() = 2.0;
        assert_eq!(c.get(0), Some(&2.0));

        // assert_eq!(c.map(|x| x * x), CartessianND::new(vec![4.0, 9.0]));

        c.map_inplace(|x| *x = *x * *x);
        assert_eq!(c, CartessianND::new(vec![4.0, 9.0]));

        let d = CartessianND::new(vec![1.0, 3.0]);
        c.zip_mut_with(&d, |x, y| {*x = *x + *y});
        assert_eq!(c, CartessianND::new(vec![5.0, 12.0]));

        assert_eq!(c.as_ref(), &[5.0, 12.0]);

        let cmut = c.as_mut();
        cmut[0] = 1.0;
        assert_eq!(c, CartessianND::new(vec![1.0, 12.0]));

        assert_eq!(c[0], 1.0);
        assert_eq!(c[1], 12.0);
        c[0] = 3.0;
        c[1] = 9.0;
        assert_eq!(c, CartessianND::new(vec![3.0, 9.0]));

        let mut iter = c.into_iter();
        assert_eq!(iter.next(), Some(&3.0));
        assert_eq!(iter.next(), Some(&9.0));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_default(){
        let  a = Cartessian2D::<i32>::default();
        assert_eq!(a, Cartessian2D::new([0, 0]));
    }
}




