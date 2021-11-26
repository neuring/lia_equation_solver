use std::{
    fmt::Display,
    ops::{AddAssign, MulAssign, RemAssign},
};

pub trait Numeric:
    Display
    + Clone
    + for<'a> AddAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + for<'a> RemAssign<&'a Self>
    + AddAssign<i64>
    + MulAssign<i64>
    + RemAssign<i64>
    + Eq
    + Ord
    + 'static
{
    /// after calling this function self contains the euclidean quotient and rhs might
    /// have changed (implementation defined)
    fn div_euc_assign(&mut self, rhs: &mut Self);

    fn gcd_assign(&mut self, rhs: &Self);
}

impl Numeric for rug::Integer {
    fn div_euc_assign(&mut self, rhs: &mut Self) {
        self.div_rem_euc_mut(rhs);
    }

    fn gcd_assign(&mut self, rhs: &Self) {
        self.gcd_mut(rhs);
    }
}

impl Numeric for i64 {
    fn gcd_assign(&mut self, rhs: &Self) {
        *self = num::integer::gcd(*self, *rhs);
    }

    fn div_euc_assign(&mut self, rhs: &mut Self) {
        *self = self.div_euclid(*rhs);
    }
}
