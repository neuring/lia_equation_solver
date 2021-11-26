use std::{
    fmt::{Display, Debug},
    ops::{AddAssign, MulAssign, RemAssign, DivAssign, SubAssign},
};

pub trait Numeric:
    Display
    + Debug
    + Clone
    + for<'a> AddAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + for<'a> RemAssign<&'a Self>
    + for<'a> DivAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + AddAssign<i64>
    + MulAssign<i64>
    + RemAssign<i64>
    + DivAssign<i64>
    + SubAssign<i64>
    + From<i64>
    + PartialOrd<i64>
    + Eq
    + PartialOrd<i64>
    + Ord
    + 'static
{
    /// after calling this function self contains the euclidean quotient and rhs might
    /// have changed (implementation defined)
    fn div_euc_assign(&mut self, rhs: &mut Self);

    /// self <- gcd(self, rhs)
    fn gcd_assign(&mut self, rhs: &Self);

    fn assign(&mut self, value: i64);

    /// self <- |value|
    fn abs_assign(&mut self, value: &Self);
}

impl Numeric for rug::Integer {
    fn div_euc_assign(&mut self, rhs: &mut Self) {
        self.div_rem_euc_mut(rhs);
    }

    fn gcd_assign(&mut self, rhs: &Self) {
        self.gcd_mut(rhs);
    }

    fn assign(&mut self, value: i64) {
        <Self as rug::Assign<i64>>::assign(self, value);
    }

    fn abs_assign(&mut self, value: &Self) {
        <Self as rug::Assign<_>>::assign(self, value.abs_ref()); 
    }
}

impl Numeric for i64 {
    fn div_euc_assign(&mut self, rhs: &mut Self) {
        *self = self.div_euclid(*rhs);
    }

    fn gcd_assign(&mut self, rhs: &Self) {
        *self = num::integer::gcd(*self, *rhs);
    }

    fn assign(&mut self, value: i64) {
        *self = value;
    }

    fn abs_assign(&mut self, value: &Self) {
        *self = value.abs();
    }
}
