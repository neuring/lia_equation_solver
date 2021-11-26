use crate::{
    math,
    numeric::Numeric,
    system::{EquationStorage, EquationViewMut, System},
};

struct ScratchData<N> {
    pub scratch_pad: Vec<N>,

    pub scratch1: N,
    pub scratch2: N,
    pub scratch3: N,
    pub scratch4: N,
}

impl<N: Numeric> ScratchData<N> {
    pub fn new(variables: usize) -> Self {
        Self {
            scratch_pad: vec![N::from(0); variables + 1],
            scratch1: N::from(0),
            scratch2: N::from(0),
            scratch3: N::from(0),
            scratch4: N::from(0),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Result {
    Sat,
    Unsat,
}

pub fn solve_equation<N: Numeric>(system: &mut System<N>) -> Result {
    let mut scratch = ScratchData::new(system.storage.variables);

    if !preprocess(system, &mut scratch) {
        return Result::Unsat;
    };

    //println!("After preprocessing\n{}", system.equations_display());

    let mut reductions_between_eliminations = 0;
    let mut eliminates_since_last_resize = 0;

    while let Some(result) =
        find_smallest_non_zero_coefficient(&system.storage, &mut scratch.scratch1)
    {
        if result.coefficient == &1 || result.coefficient == &-1 {
            eliminate_equation(
                system,
                result.equation_idx,
                result.coefficient_idx,
                &mut scratch,
            );
            println!(
                "eliminated: {} (reductions inbetween: {})",
                system.killed_variables, reductions_between_eliminations,
            );
            reductions_between_eliminations = 0;
            //println!("After elimination\n{}", system.equations_display());
            eliminates_since_last_resize += 1;

            if eliminates_since_last_resize > 20 && system.storage.variables > 4 {
                //println!("resizing...");
                //dbg!(&system.alive_terms);
                //dbg!(&system.varmap);
                //dbg!(&system.storage);
                system.resize();
                eliminates_since_last_resize = 0;
                //dbg!(&system.varmap);
                //dbg!(&system.storage);
            }
        } else {
            assert_ne!(result.coefficient, &0);
            reduce_coefficients(
                system,
                result.equation_idx,
                result.coefficient_idx,
                &mut scratch,
            );

            //println!("After reduce\n{}", system.equations_display());

            reductions_between_eliminations += 1;
        }
        if !preprocess(system, &mut scratch) {
            return Result::Unsat;
        };
        //println!("After simplification\n{}", system.equations_display());
    }

    if find_any_contradictions(system.get_storage()) {
        Result::Unsat
    } else {
        Result::Sat
    }
}

fn preprocess<N: Numeric>(
    system: &mut System<N>,
    scratch: &mut ScratchData<N>,
) -> bool {
    let storage = &mut system.storage;
    for mut equation in storage.iter_equations_mut().filter(|eq| !eq.is_empty()) {
        let gcd = &mut scratch.scratch1;
        math::gcd(gcd, equation.get_coefficient_slice());

        let rem = &mut scratch.scratch2;
        rem.clone_from(equation.get_result());
        *rem %= &*gcd;

        if rem.cmp_zero().is_ne() {
            return false;
        }

        if *gcd > 1 {
            // TODO: does this improve performance
            for coefficient in equation.iter_coefficients() {
                *coefficient /= &*gcd;
            }
        }

        *equation.get_result() /= &*gcd;
    }

    true
}

#[derive(Debug, Clone, Copy)]
struct SearchResult<'a, N> {
    equation_idx: usize,
    coefficient_idx: usize,
    coefficient: &'a N,
}

fn find_smallest_non_zero_coefficient<'a, N: Numeric>(
    system: &EquationStorage<N>,
    scratch1: &'a mut N,
) -> Option<SearchResult<'a, N>> {
    let mut found_min = false;

    let current_min = scratch1;
    let mut current_coefficient_idx = 0;
    let mut current_equation_idx = 0;

    let mut one_coef_counter = 0;

    for (equation_idx, equation) in system.iter_equations().enumerate() {
        for (coefficient_idx, coefficient) in
            equation.iter_coefficients().enumerate()
        {
            if coefficient == &1 || coefficient == &-1 {
                one_coef_counter += 1;
            }

            let is_new_minimum = coefficient.cmp_zero().is_ne()
                && (!found_min || coefficient.abs_compare(&current_min).is_lt());

            if is_new_minimum {
                found_min = true;
                current_min.clone_from(&coefficient);
                current_coefficient_idx = coefficient_idx;
                current_equation_idx = equation_idx;
            }
        }
    }

    //if one_coef_counter > 1 {
    //    println!("Had elimination choice {}", one_coef_counter);
    //}

    found_min.then(move || SearchResult {
        equation_idx: current_equation_idx,
        coefficient_idx: current_coefficient_idx,
        coefficient: &*current_min,
    })
}

fn eliminate_equation<N: Numeric>(
    system: &mut System<N>,
    eliminated_equation_idx: usize,
    eliminated_coefficient_idx: usize,
    scratch: &mut ScratchData<N>,
) {
    let storage = &mut system.storage;
    let varmap = &system.varmap;

    let mut eliminated_equation = EquationViewMut {
        data: &mut scratch.scratch_pad,
    };
    eliminated_equation.copy_into(storage.get_equation(eliminated_equation_idx));
    let eliminated_equation = eliminated_equation.into_ref();

    let eliminated_coefficient =
        eliminated_equation.get_coefficient(eliminated_coefficient_idx);

    let sign = if eliminated_coefficient.cmp_zero().is_gt() {
        1
    } else {
        -1
    };

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
                .map(|(i, c)| {
                    let mut result = N::from(-sign);
                    result *= &*c;
                    (varmap[i], result)
                })
                .collect();

            let mut constant = N::from(sign);
            constant *= &*equation.get_result();
            system.reconstruction.add(var, terms, constant);

            equation.clear();
        } else {
            let target_coefficient_factor = &mut scratch.scratch2;
            target_coefficient_factor
                .clone_from(equation.get_coefficient(eliminated_coefficient_idx));

            for (target_coefficient, coefficient) in equation
                .iter_coefficients()
                .zip(eliminated_equation.iter_coefficients())
            {
                let s = &mut scratch.scratch1;
                s.clone_from(&*target_coefficient_factor);
                *s *= sign;
                *s *= coefficient;

                *target_coefficient -= &*s;
            }

            equation
                .get_coefficient(eliminated_coefficient_idx)
                .assign(0);

            let s = &mut scratch.scratch1;
            s.clone_from(&*target_coefficient_factor);
            *s *= sign;
            *s *= eliminated_equation.get_result();

            *equation.get_result() -= &*s;
        }
    }

    system.kill_variable(eliminated_coefficient_idx);
}

fn reduce_coefficients<N: Numeric>(
    system: &mut System<N>,
    equation_idx: usize,
    coefficient_idx: usize,
    scratch: &mut ScratchData<N>,
) {
    let storage = &mut system.storage;

    let mut equation = storage.get_equation_mut(equation_idx);

    let coefficient = &mut scratch.scratch4;
    coefficient.clone_from(equation.get_coefficient(coefficient_idx));

    // println!("reducing coefficient {}", coefficient);

    let original_coefficient_idx = coefficient_idx;

    if coefficient.cmp_zero().is_lt() {
        equation.iter_coefficients().for_each(|c| *c *= -1);
        *equation.get_result() *= -1;
        *coefficient *= -1;
    }

    let m = &mut scratch.scratch1;
    m.clone_from(&*coefficient);
    *m += 1;

    for (coefficient_idx, coefficient) in equation.iter_coefficients().enumerate() {
        if coefficient_idx == original_coefficient_idx {
            *coefficient *= -1;

            let minus_m = &mut scratch.scratch2;
            minus_m.clone_from(m);
            *minus_m *= -1;
            scratch.scratch_pad[coefficient_idx].clone_from(&*minus_m);
        } else {
            let rounded_div = &mut scratch.scratch2;
            let sm = &mut scratch.scratch3;
            math::special_mod(sm, rounded_div, &*coefficient, &m);

            coefficient.clone_from(&rounded_div);
            *coefficient += &*sm;
            scratch.scratch_pad[coefficient_idx].clone_from(sm);
        }
    }

    let result = equation.get_result();
    let rounded_div = &mut scratch.scratch2;
    let sm = &mut scratch.scratch3;
    math::special_mod(sm, rounded_div, &*result, &*m);

    result.clone_from(&rounded_div);
    *result += &*sm;

    let equation_sm = sm;

    for (e_idx, mut equation) in storage
        .iter_equations_mut()
        .enumerate()
        .filter(|(_, eq)| !eq.is_empty())
    {
        if equation_idx != e_idx {
            let coefficient_factor = &mut scratch.scratch4;
            coefficient_factor.clone_from(equation.get_coefficient(coefficient_idx));

            for (c_idx, coefficient) in equation.iter_coefficients().enumerate() {
                if c_idx == coefficient_idx {
                    coefficient.clone_from(&*m);
                    *coefficient *= -1;
                    *coefficient *= &*coefficient_factor;
                } else {
                    let summand = &mut scratch.scratch2;
                    summand.clone_from(&coefficient_factor);
                    *summand *= &scratch.scratch_pad[c_idx];
                    *coefficient += &*summand;
                }
            }

            let summand = &mut scratch.scratch2;
            summand.clone_from(&coefficient_factor);
            *summand *= &*equation_sm;
            *equation.get_result() += &*summand;
        }
    }

    let old_var = system.varmap[original_coefficient_idx];
    let new_var = system.new_variable();
    system.map_variable(original_coefficient_idx, new_var);

    //println!(
    //    "{} = {}",
    //    util::fmt_variable(old_var, system.storage.variables),
    //    itertools::join(
    //        scratch
    //            .scratch_pad
    //            .iter()
    //            .enumerate()
    //            .filter(|(_, c)| *c != &0)
    //            .map(|(i, c)| (system.varmap[i], c))
    //            .map(|(x, c)| format!(
    //                "{}{}",
    //                c,
    //                util::fmt_variable(x, system.storage.variables)
    //            ))
    //            .chain(std::iter::once(format!("{}", {
    //                let mut x = equation_sm.clone();
    //                x *= -1;
    //                x
    //            }))),
    //        " + "
    //    )
    //);

    let terms = scratch
        .scratch_pad
        .iter()
        .take(system.storage.variables)
        .enumerate()
        .filter(|(_, c)| c.cmp_zero().is_ne())
        .map(|(i, c)| (system.varmap[i], c.clone()))
        .collect();
    let mut constant = equation_sm.clone();
    constant *= -1;
    system.reconstruction.add(old_var, terms, constant);
}

fn find_any_contradictions<N: Numeric>(storage: &EquationStorage<N>) -> bool {
    storage.iter_equations().any(|eq| !eq.is_empty())
}
