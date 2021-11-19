
pub fn subscript(mut i: u32) -> String {
    let mut subscript_digits = Vec::new();

    if i == 0  {
        subscript_digits.push('₀');
    }

    while i != 0 {
        let digit = i % 10;
        i /= 10;

        subscript_digits.push(char::from_u32(('₀' as u32) + digit).unwrap());
    }

    subscript_digits.reverse();
    subscript_digits.into_iter().collect()
}
