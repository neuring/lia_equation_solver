use std::{
    fmt::{Debug, Display},
    ops::{AddAssign, DivAssign, MulAssign, RemAssign, SubAssign},
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

    fn abs_compare(&self, value: &Self) -> std::cmp::Ordering;
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

    fn abs_compare(&self, value: &Self) -> std::cmp::Ordering {
        self.cmp_abs(value)
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

    fn abs_compare(&self, value: &Self) -> std::cmp::Ordering {
        self.abs().cmp(&value.abs())
    }
}
