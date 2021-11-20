use crate::system::VariableIndex;

pub fn fmt_variable(i: VariableIndex, total_variables: usize) -> String {
    if i.0 < total_variables {
        format!("x{}", subscript(i.0 as _))
    } else {
        format!("y{}", subscript(i.0 as _))
    }
}

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
