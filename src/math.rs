use num::Integer;

pub fn gcd(values: &[i32]) -> i32 {
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

pub fn special_mod(a: i32, b: i32) -> i32 {
    a - b * ((a as f64) / (b as f64) + 0.5) as i32
}

/// calculates ⌊a / b + 1/2 ⌋
pub fn rounded_divisor(a: i32, b: i32) -> i32 {
    ((a as f64) / (b as f64) + 0.5) as i32
}
