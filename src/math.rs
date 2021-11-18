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
