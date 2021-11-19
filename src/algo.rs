use crate::{math, system::System};

pub enum Result {
    Sat,
    Unsat,
}

pub fn solve_equation(system: &mut System) -> Result {
    if !preprocess(system) {
        return Result::Unsat;
    };

    while let Some(result) = find_smallest_non_zero_coefficient(system) {
        if result.coefficient.abs() == 1 {
            eliminate_equation(system, result.equation_idx, result.coefficient_idx);
        } else {
            assert!(result.coefficient != 0);
            reduce_coefficients(system, result.equation_idx, result.coefficient_idx);
        }
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

#[derive(Debug, Clone, Copy)]
struct SearchResult {
    equation_idx: usize,
    coefficient_idx: usize,
    coefficient: i32,
}

fn find_smallest_non_zero_coefficient(system: &System) -> Option<SearchResult> {
    let mut min: Option<SearchResult> = None;

    for (equation_idx, equation) in system.iter_equations().enumerate() {
        for (coefficient_idx, coefficient) in
            equation.iter_coefficients().enumerate()
        {
            let new_minimum = coefficient.abs() > 0
                && min
                    .map(|min| min.coefficient.abs() > coefficient.abs())
                    .unwrap_or(true);

            if new_minimum {
                min = Some(SearchResult {
                    equation_idx,
                    coefficient_idx,
                    coefficient,
                });
            }
        }
    }

    min
}

fn eliminate_equation(
    system: &mut System,
    eliminated_equation_idx: usize,
    eliminated_coefficient_idx: usize,
) {
    let eliminated_equation =
        system.get_equation_mut(eliminated_equation_idx).to_owned();

    let eliminated_coefficient =
        eliminated_equation.get_coefficient(eliminated_coefficient_idx);

    let sign = if eliminated_coefficient > 0 { 1 } else { -1 };

    for (equation_idx, mut equation) in system.iter_equations_mut().enumerate() {
        if equation_idx == eliminated_equation_idx {
            equation.clear();
        } else {
            for (target_coefficient, coefficient) in equation
                .iter_coefficients()
                .zip(eliminated_equation.iter_coefficients())
            {
                *target_coefficient -= sign * coefficient;
            }

            *equation.get_coefficient(eliminated_coefficient_idx) = 0;

            *equation.get_result() -= sign * eliminated_equation.get_result();
        }
    }
}

fn reduce_coefficients(
    system: &mut System,
    equation_idx: usize,
    coefficient_idx: usize,
) {
    let mut equation = system.get_equation_mut(equation_idx);

    let mut coefficient = *equation.get_coefficient(coefficient_idx);

    let original_coefficient_idx = coefficient_idx;

    if coefficient < 0 {
        equation.iter_coefficients().for_each(|c| *c *= -1);
        coefficient *= -1;
    }

    let m = coefficient + 1;

    for (coefficient_idx, coefficient) in equation.iter_coefficients().enumerate() {
        if coefficient_idx == original_coefficient_idx {
            *coefficient = -1;
        } else {
            *coefficient = math::rounded_divisor(*coefficient, m)
                + math::special_mod(*coefficient, m);
        }
    }

    let result = equation.get_result();
    *result = math::rounded_divisor(*result, m) + math::special_mod(*result, m);
}
