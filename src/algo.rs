use crate::{math, system::System};

pub enum Result {
    Sat,
    Unsat,
}

pub fn solve_equation(system: &mut System) -> Result {

    if !preprocess(system) { return Result::Unsat };

    while let Some((equation_idx, coefficient_idx)) = can_eliminate_equation(system) {
        eliminate_equation(system, equation_idx, coefficient_idx);
    }

    Result::Sat
}

fn preprocess(system: &mut System) -> bool {
    for mut equation in system.iter_equations_mut() {

        let gcd = math::gcd(equation.get_coefficient_slice());

        if *equation.get_result() % gcd != 0 {
            return false;
        }

        for coefficient in equation.iter_coefficients() {
            *coefficient /= gcd;
        }

        *equation.get_result() /= gcd;
    }

    true
}

fn can_eliminate_equation(system: &System) -> Option<(usize, usize)> {
    for (equation_idx, equation) in system.iter_equations().enumerate() {
        for (coefficient_idx, coefficient) in
            equation.iter_coefficients().enumerate()
        {
            if coefficient.abs() == 1 {
                return Some((equation_idx, coefficient_idx));
            }
        }
    }

    None
}

fn eliminate_equation(
    system: &mut System,
    eliminated_equation_idx: usize,
    eliminated_coefficient_idx: usize,
) {
    let eliminated_equation = system.get_equation_mut(eliminated_equation_idx).to_owned();

    let eliminated_coefficient = eliminated_equation.get_coefficient(eliminated_coefficient_idx);

    let sign = if eliminated_coefficient > 0 { 1 } else { -1 };

    for (equation_idx, mut equation) in system.iter_equations_mut().enumerate() {

        if equation_idx == eliminated_equation_idx {
            equation.clear();
        } else {

            for (target_coefficient, coefficient) in equation.iter_coefficients().zip(eliminated_equation.iter_coefficients()) {
                *target_coefficient -= sign * coefficient;
            }

            *equation.get_coefficient(eliminated_coefficient_idx) = 0;

            *equation.get_result() -= sign * eliminated_equation.get_result();
        }
    }
}
