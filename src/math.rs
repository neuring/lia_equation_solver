use num::Integer;

pub fn gcd(values: &[i64]) -> i64 {
    if values.len() == 1 {
        values[0]
    } else {
        let half = values.len() / 2;

        let left = &values[..half];
        let right = &values[half..];

        let left_gcd = gcd(left);
        let right_gcd = gcd(right);

        left_gcd.gcd(&right_gcd)
    }
}

pub fn special_mod(a: i64, b: i64) -> i64 {
    a - b * f64::floor((a as f64) / (b as f64) + 0.5) as i64
}

/// calculates ⌊a / b + 1/2 ⌋
pub fn rounded_divisor(a: i64, b: i64) -> i64 {
    f64::floor((a as f64) / (b as f64) + 0.5) as i64
}
