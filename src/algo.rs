use crate::{
    math,
    system::{Equation, EquationStorage, System},
    util,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Result {
    Sat,
    Unsat,
}

fn free_space(eliminated: usize, variables: usize, equations: usize) -> usize {
    eliminated * (variables * 1) + (equations - eliminated) * eliminated
}

pub fn solve_equation(system: &mut System) -> Result {
    if !preprocess(system.get_storage_mut()) {
        return Result::Unsat;
    };

    println!("After preprocessing\n{}", system.equations_display());

    let mut eliminations = 0;
    let mut reductions_between_eliminations = 0;
    while let Some(result) = find_smallest_non_zero_coefficient(system.get_storage())
    {
        if result.coefficient.abs() == 1 {
            eliminate_equation(system, result.equation_idx, result.coefficient_idx);
            eliminations += 1;
            println!(
                "eliminated: {} (reductions inbetween: {}, estimate free space: {}, actual: {})",
                eliminations,
                reductions_between_eliminations,
                free_space(
                    eliminations,
                    system.storage.variables,
                    system.storage.equations
                ),
                system.storage.data.iter().filter(|&&i| i == 0).count()
            );
            reductions_between_eliminations = 0;
            println!("After elimination\n{}", system.equations_display());
        } else {
            assert!(result.coefficient != 0);
            reduce_coefficients(system, result.equation_idx, result.coefficient_idx);
            println!("After reduce\n{}", system.equations_display());
            reductions_between_eliminations += 1;
        }
        if !preprocess(system.get_storage_mut()) {
            return Result::Unsat;
        };
        println!("After simplification\n{}", system.equations_display());
    }

    if find_any_contradictions(system.get_storage()) {
        Result::Unsat
    } else {
        Result::Sat
    }
}

fn preprocess(system: &mut EquationStorage) -> bool {
    for mut equation in system.iter_equations_mut().filter(|eq| !eq.is_empty()) {
        let mut gcd = 0;
        math::gcd(&mut gcd, equation.get_coefficient_slice());

        if *equation.get_result() % gcd != 0 {
            return false;
        }

        if gcd > 1 {
            // TODO: does this improve performance
            for coefficient in equation.iter_coefficients() {
                *coefficient /= gcd;
            }
        }

        *equation.get_result() /= gcd;
    }

    true
}

#[derive(Debug, Clone, Copy)]
struct SearchResult {
    equation_idx: usize,
    coefficient_idx: usize,
    coefficient: i64,
}

fn find_smallest_non_zero_coefficient(
    system: &EquationStorage,
) -> Option<SearchResult> {
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
    let storage = &mut system.storage;
    let varmap = &system.varmap;
    // TODO: copy in scratch_pad instead of allocating
    let eliminated_equation =
        storage.get_equation_mut(eliminated_equation_idx).to_owned();

    let eliminated_coefficient =
        eliminated_equation.get_coefficient(eliminated_coefficient_idx);

    let sign = if eliminated_coefficient > 0 { 1 } else { -1 };

    for (equation_idx, mut equation) in storage
        .iter_equations_mut()
        .enumerate()
        .filter(|(_, eq)| !eq.is_empty())
    {
        if equation_idx == eliminated_equation_idx {
            // Add eliminated equation to reconstruction
            let var = system.varmap[eliminated_coefficient_idx];
            let terms = equation
                .iter_coefficients()
                .enumerate()
                .filter(|(i, c)| **c != 0 && *i != eliminated_coefficient_idx)
                .map(|(i, c)| (varmap[i], -sign * *c))
                .collect();
            system
                .reconstruction
                .add(var, terms, sign * *equation.get_result());

            equation.clear();
        } else {
            let target_coefficient_factor =
                *equation.get_coefficient(eliminated_coefficient_idx);

            for (target_coefficient, coefficient) in equation
                .iter_coefficients()
                .zip(eliminated_equation.iter_coefficients())
            {
                *target_coefficient -=
                    sign * target_coefficient_factor * coefficient;
            }

            *equation.get_coefficient(eliminated_coefficient_idx) = 0;

            *equation.get_result() -=
                sign * target_coefficient_factor * eliminated_equation.get_result();
        }
    }
}

fn reduce_coefficients(
    system: &mut System,
    equation_idx: usize,
    coefficient_idx: usize,
) {
    let storage = &mut system.storage;

    let mut equation = storage.get_equation_mut(equation_idx);

    let mut coefficient = *equation.get_coefficient(coefficient_idx);
    // println!("reducing coefficient {}", coefficient);

    let original_coefficient_idx = coefficient_idx;

    if coefficient < 0 {
        equation.iter_coefficients().for_each(|c| *c *= -1);
        *equation.get_result() *= -1;
        coefficient *= -1;
    }

    let m = coefficient + 1;

    for (coefficient_idx, coefficient) in equation.iter_coefficients().enumerate() {
        if coefficient_idx == original_coefficient_idx {
            *coefficient *= -1;

            system.scratch_pad[coefficient_idx] = -m;
        } else {
            let mut rounded_div = 0;
            let mut sm = 0;
            math::special_mod(&mut sm, &mut rounded_div, coefficient, &m);

            *coefficient = rounded_div + sm;
            system.scratch_pad[coefficient_idx] = sm;
        }
    }

    let result = equation.get_result();
    let mut rounded_div = 0;
    let mut sm = 0;
    math::special_mod(&mut sm, &mut rounded_div, result, &m);

    *result = rounded_div  + sm;

    let equation_sm = sm;

    for (e_idx, mut equation) in storage.iter_equations_mut().enumerate() {
        if equation_idx != e_idx {
            let coefficient_factor = *equation.get_coefficient(coefficient_idx);

            for (c_idx, coefficient) in equation.iter_coefficients().enumerate() {
                if c_idx == coefficient_idx {
                    *coefficient = -m * coefficient_factor;
                } else {
                    *coefficient += coefficient_factor * system.scratch_pad[c_idx];
                }
            }

            *equation.get_result() += coefficient_factor * equation_sm;
        }
    }

    let old_var = system.varmap[original_coefficient_idx];
    let new_var = system.new_variable();
    system.map_variable(original_coefficient_idx, new_var);

    /*println!(
        "{} = {}",
        util::fmt_variable(old_var, system.storage.variables),
        itertools::join(
            system
                .scratch_pad
                .iter()
                .enumerate()
                .filter(|&(_, &c)| c != 0)
                .map(|(i, c)| (system.varmap[i], *c))
                .map(|(x, c)| format!(
                    "{}{}",
                    c,
                    util::fmt_variable(x, system.storage.variables)
                ))
                .chain(std::iter::once(format!("{}", -equation_sm))),
            " + "
        )
    );*/

    let terms = system
        .scratch_pad
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c != 0)
        .map(|(i, c)| (system.varmap[i], *c))
        .collect();
    system.reconstruction.add(old_var, terms, -equation_sm);
}

fn find_any_contradictions(storage: &EquationStorage) -> bool {
    storage.iter_equations().any(|eq| !eq.is_empty())
}
