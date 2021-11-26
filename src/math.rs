use crate::numeric::Numeric;

pub fn gcd<N: Numeric>(result: &mut N, values: &[N]) {
    result.clone_from(&values[0]);

    for val in values.iter().skip(1) {
        result.gcd_assign(val);
    }
}

pub fn special_mod<N: Numeric>(
    result: &mut N,
    rounded_div: &mut N,
    a: &N,
    b: &N,
) {
    rounded_divisor(rounded_div, a, b, result);

    result.clone_from(rounded_div);
    *result *= b;
    *result *= -1;
    *result += a;
}

/// calculates ⌊a / b + 1/2 ⌋, but without using floats
/// ⌊a / b + 1/2 ⌋ = ⌊(a + 2 * b) / (2 * b)⌋
pub fn rounded_divisor<N: Numeric>(
    result: &mut N,
    a: &N,
    b: &N,
    scratch: &mut N,
) {
    result.clone_from(a);
    *result *= 2;
    *result += b;

    scratch.clone_from(b);
    *scratch *= 2;

    result.div_euc_assign(scratch);
}
