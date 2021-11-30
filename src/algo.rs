use crate::{
    math,
    numeric::Numeric,
    system::{EquationStorage, EquationView, EquationViewMut, System},
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

    //let mut debug_out = std::io::BufWriter::new(
    //    std::fs::File::create("simplifications_per_iteration").unwrap(),
    //);

    if !preprocess(system, &mut scratch) {
        return Result::Unsat;
    };

    //println!("After preprocessing\n{}", system.equations_display());

    let mut reductions_between_eliminations = 0;
    let mut eliminates_since_last_resize = 0;

    while let Some(result) =
        find_smallest_non_zero_coefficient(&system.storage, &mut scratch)
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
                system.resize();
                eliminates_since_last_resize = 0;
            }

            //if system.killed_variables % 1 == 0 {
            //    system
            //        .storage
            //        .print_value_stats(system.killed_variables, &mut debug_out);
            //}
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
        if !simplify(system, &mut scratch) {
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

        if *gcd > 1 {
            let rem = &mut scratch.scratch2;
            rem.clone_from(equation.get_result());
            rem.negate();
            *rem %= &*gcd;

            if rem.cmp_zero().is_ne() {
                return false;
            }

            for coefficient in equation.iter_coefficients() {
                *coefficient /= &*gcd;
            }

            *equation.get_result() /= &*gcd;
        }
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
    scratch: &'a mut ScratchData<N>,
) -> Option<SearchResult<'a, N>> {
    let mut found_min = false;

    let current_min = &mut scratch.scratch1;
    let mut min_equation_idx = 0;
    let mut min_coefficient_idx = 0;

    for (equation_idx, equation) in system
        .iter_equations()
        .enumerate()
        .filter(|(_, eq)| !eq.is_empty())
    {
        for (coefficient_idx, coefficient) in equation
            .iter_coefficients()
            .enumerate()
            .filter(|(_, coef)| coef.cmp_zero().is_ne())
        {
            if !found_min || coefficient.abs_compare(&current_min).is_lt() {
                found_min = true;
                current_min.clone_from(&coefficient);
                min_equation_idx = equation_idx;
                min_coefficient_idx = coefficient_idx;
            }
        }
    }

    found_min.then(move || SearchResult {
        equation_idx: min_equation_idx,
        coefficient_idx: min_coefficient_idx,
        coefficient: current_min,
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

    let storage_variables = storage.variables;

    let mut eliminated_equation = storage.get_equation_mut(eliminated_equation_idx);
    let eliminated_coefficient =
        eliminated_equation.get_coefficient(eliminated_coefficient_idx);

    // Normalise target coefficient to be 1 (if it is -1)
    if eliminated_coefficient.cmp_zero().is_lt() {
        eliminated_equation
            .data
            .iter_mut()
            .for_each(|value| value.negate());
    }

    let mut eliminated_equation_copy = EquationViewMut {
        data: &mut scratch.scratch_pad[..storage_variables + 1],
    };
    eliminated_equation_copy.copy_into(eliminated_equation.as_ref());

    // Remove equation from storage
    eliminated_equation.clear();

    // Make equation to term which can be substituted into all other equations.
    eliminated_equation_copy
        .get_coefficient(eliminated_coefficient_idx)
        .assign(0);

    eliminated_equation_copy
        .data
        .iter_mut()
        .for_each(|v| v.negate());

    // seal term to prevent accidental mutation
    let eliminated_equation = eliminated_equation_copy.into_ref();

    // Add eliminated equation to reconstruction
    let terms = eliminated_equation
        .iter_coefficients()
        .enumerate()
        .filter(|&(_, c)| c.cmp_zero().is_ne())
        .map(|(i, c)| (varmap[i], c.clone()))
        .collect();

    let constant = eliminated_equation.get_result().clone();

    let var = system.varmap[eliminated_coefficient_idx];
    system.reconstruction.add(var, terms, constant);

    for equation in storage.iter_equations_mut().filter(|eq| !eq.is_empty()) {
        substitute_variable_with_term(
            equation,
            eliminated_coefficient_idx,
            eliminated_equation,
            &mut scratch.scratch3,
            &mut scratch.scratch4,
        );
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

    let storage_variables = storage.variables;

    let mut equation = storage.get_equation_mut(equation_idx);

    let coefficient = &mut scratch.scratch4;
    coefficient.clone_from(equation.get_coefficient(coefficient_idx));

    // println!("reducing coefficient {}", coefficient);

    let original_coefficient_idx = coefficient_idx;

    if coefficient.cmp_zero().is_lt() {
        equation.data.iter_mut().for_each(|c| c.negate());
        coefficient.negate();
    }

    let m = &mut scratch.scratch1;
    m.clone_from(&*coefficient);
    *m += 1;

    let substitution_term = EquationViewMut {
        data: &mut scratch.scratch_pad[..storage_variables + 1],
    };

    for (coefficient, coefficient_idx) in equation.data.iter_mut().zip(0..) {
        // + 1, because this loop starts with the result element as zero, not the first coefficient.
        if coefficient_idx == original_coefficient_idx + 1 {
            assert!(coefficient_idx > 0);

            coefficient.negate();

            let minus_m = &mut scratch.scratch2;
            minus_m.clone_from(m);
            minus_m.negate();

            substitution_term.data[coefficient_idx].clone_from(&*minus_m);
        } else {
            let rounded_div = &mut scratch.scratch2;
            let sm = &mut scratch.scratch3;
            math::special_mod(sm, rounded_div, &*coefficient, &m);

            coefficient.clone_from(&rounded_div);
            *coefficient += &*sm;
            substitution_term.data[coefficient_idx].clone_from(sm);
        }
    }

    let substitution_term = substitution_term.into_ref();

    for equation in storage
        .iter_equations_mut()
        .enumerate()
        .filter(|(eq_idx, eq)| !eq.is_empty() && *eq_idx != equation_idx)
        .map(|(_, eq)| eq)
    {
        substitute_variable_with_term(
            equation,
            original_coefficient_idx,
            substitution_term,
            &mut scratch.scratch2,
            &mut scratch.scratch3,
        )
    }

    let old_var = system.varmap[original_coefficient_idx];
    let new_var = system.new_variable();
    system.map_variable(original_coefficient_idx, new_var);

    // Add substitution to reconstruction
    let terms = substitution_term
        .iter_coefficients()
        .enumerate()
        .filter(|(_, c)| c.cmp_zero().is_ne())
        .map(|(i, c)| (system.varmap[i], c.clone()))
        .collect();
    let constant = substitution_term.get_result().clone();
    system.reconstruction.add(old_var, terms, constant);
}

fn simplify<N: Numeric>(
    system: &mut System<N>,
    scratch: &mut ScratchData<N>,
) -> bool {
    let storage = &mut system.storage;

    'outer: for mut equation in storage.iter_equations_mut() {
        let mut non_zero_coef = None;

        for (coefficient_idx, coefficient) in
            equation.iter_coefficients().enumerate()
        {
            if coefficient.cmp_zero().is_ne() {
                if non_zero_coef.is_some() {
                    continue 'outer;
                } else {
                    non_zero_coef = Some(coefficient_idx);
                }
            }
        }

        if let Some(idx) = non_zero_coef {
            let coef = &mut scratch.scratch1;
            coef.assign(1);

            std::mem::swap(coef, equation.get_coefficient(idx));

            let modulo = &mut scratch.scratch2;
            modulo.clone_from(equation.get_result());
            modulo.negate();
            *modulo %= &*coef;

            if modulo.cmp_zero().is_ne() {
                return false;
            }

            *equation.get_result() /= &*coef;
        }
    }

    true
}

fn find_any_contradictions<N: Numeric>(storage: &EquationStorage<N>) -> bool {
    storage.iter_equations().any(|eq| !eq.is_empty())
}

// Inserts `term` at coefficient index `target_idx`.
// The term might have a new variable at `target_idx` which then replaces
// the variable originally in `self`
fn substitute_variable_with_term<N: Numeric>(
    mut target: EquationViewMut<'_, N>,
    target_idx: usize,
    term: EquationView<'_, N>,
    scratch1: &mut N,
    scratch2: &mut N,
) {
    let factor = scratch1;

    let target_coef = target.get_coefficient(target_idx);

    factor.clone_from(target_coef);

    target_coef.assign(0);

    let result = scratch2;
    for (target_coef, term_coef) in target.data.iter_mut().zip(term.data.iter()) {
        result.clone_from(factor);
        *result *= term_coef;
        *target_coef += &*result;
    }
}
