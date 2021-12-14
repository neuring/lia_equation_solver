use crate::{
    math,
    numeric::Numeric,
    system::{Equation, EquationStorage, EquationView, EquationViewMut, System},
    SolverConfig,
};

use rayon::iter::{IndexedParallelIterator, ParallelIterator};

struct ScratchData<N> {
    pub scratch_pad: Vec<N>,

    pub scratch1: N,
    pub scratch2: N,
    pub scratch3: N,
    pub scratch4: N,
    pub scratch5: N,
}

impl<N: Numeric> ScratchData<N> {
    pub fn new(variables: usize) -> Self {
        Self {
            scratch_pad: vec![N::from(0); variables + 1],
            scratch1: N::from(0),
            scratch2: N::from(0),
            scratch3: N::from(0),
            scratch4: N::from(0),
            scratch5: N::from(0),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Result {
    Sat,
    Unsat,
}

pub fn solve_equation<N: Numeric>(
    system: &mut System<N>,
    config: &SolverConfig,
) -> Result {
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
                config,
            );

            if !config.silent {
                println!(
                    "eliminated: {} (reductions inbetween: {})",
                    system.killed_variables, reductions_between_eliminations,
                );
            }
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
            reductions_between_eliminations = reduce_coefficients(
                system,
                result.equation_idx,
                &mut scratch,
                config,
            );

            //println!("After reduce\n{}", system.equations_display());
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
            if !found_min || coefficient.abs_compare(current_min).is_lt() {
                found_min = true;
                current_min.clone_from(coefficient);
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
    config: &SolverConfig,
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
    eliminated_equation_copy.copy_from(eliminated_equation.as_ref());

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

    if !config.no_solution {
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
    }

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
    scratch: &mut ScratchData<N>,
    config: &SolverConfig,
) -> u32 {
    let storage = &mut system.storage;

    let storage_variables = storage.variables;

    let mut new_varmap = system.varmap.clone();

    let mut substitutions: Vec<Option<Equation<N>>> = vec![None; storage_variables];

    let mut equation = storage.get_equation_mut(equation_idx);

    //println!("reducing equation {}", equation.display(&new_varmap));

    let mut reductions = 0;
    loop {
        //println!("reduce");
        if equation
            .iter_coefficients()
            .filter(|c| c.cmp_zero().is_ne())
            .count()
            <= 1
        {
            break;
        }

        let cur_min_coef = &mut scratch.scratch1;
        let (cur_min_coef_idx, min_coef) = equation
            .iter_coefficients()
            .enumerate()
            .filter(|(_, v)| v.cmp_zero().is_ne())
            .min_by(|(_, a), (_, b)| a.abs_compare(b))
            .unwrap();

        //println!("min_coef = {}", min_coef);

        if min_coef == &1 || min_coef == &-1 {
            break;
        }

        cur_min_coef.clone_from(&*min_coef);

        if cur_min_coef.cmp_zero().is_lt() {
            equation.data.iter_mut().for_each(|c| c.negate());
            cur_min_coef.negate();
        }

        let m = &mut scratch.scratch2;
        m.clone_from(&*cur_min_coef);
        *m += 1;

        let minus_m = &mut scratch.scratch3;
        minus_m.clone_from(m);
        minus_m.negate();

        let mut substitution_term = EquationViewMut {
            data: &mut scratch.scratch_pad[..storage_variables + 1],
        };

        for (coefficient, coefficient_idx) in equation.data.iter_mut().zip(0..) {
            // + 1, because this loop starts with the result element as zero, not the first coefficient.
            if coefficient_idx == cur_min_coef_idx + 1 {
                assert!(coefficient_idx > 0);

                coefficient.negate();

                substitution_term.data[coefficient_idx].clone_from(&*minus_m);
            } else {
                let rounded_div = &mut scratch.scratch4;
                let sm = &mut scratch.scratch5;
                math::special_mod(sm, rounded_div, &*coefficient, m);

                coefficient.clone_from(rounded_div);
                *coefficient += &*sm;
                substitution_term.data[coefficient_idx].clone_from(sm);
            }
        }

        for (idx, subst) in substitutions.iter_mut().enumerate() {
            match subst {
                Some(subst) => substitute_variable_with_term(
                    subst.as_mut(),
                    cur_min_coef_idx,
                    substitution_term.as_ref(),
                    &mut scratch.scratch4,
                    &mut scratch.scratch5,
                ),
                None if idx == cur_min_coef_idx => {
                    *subst = Some(substitution_term.to_owned())
                }
                _ => {}
            }
        }

        let old_var = new_varmap[cur_min_coef_idx];
        let new_var = system.var_generator.new_variable();
        new_varmap[cur_min_coef_idx] = new_var;

        //println!("equation after reduction {}", equation.display(&new_varmap));

        //println!(
        //    "new substitution {} = {}",
        //    crate::util::fmt_variable(old_var, storage_variables),
        //    substitution_term.display(&new_varmap)
        //);

        //println!("collected substitutions:");
        //for (idx, s) in substitutions.iter().enumerate() {
        //    if let Some(s) = s {
        //        println!(
        //            "{} = {}",
        //            crate::util::fmt_variable(system.varmap[idx], storage_variables),
        //            s.display(&new_varmap)
        //        );
        //    }
        //}

        // Add substitution to reconstruction
        if !config.no_solution {
            let terms = substitution_term
                .iter_coefficients()
                .enumerate()
                .filter(|(_, c)| c.cmp_zero().is_ne())
                .map(|(i, c)| (new_varmap[i], c.clone()))
                .collect();
            let constant = substitution_term.get_result().clone();
            system.reconstruction.add(old_var, terms, constant);
        }

        reductions += 1;
    }

    //println!("apply");

    storage
        .par_iter_equations_mut()
        .enumerate()
        .filter(|(eq_idx, _)| *eq_idx != equation_idx)
        .for_each(|(_, equation)| {
            apply_all_substitutions(equation, &substitutions, storage_variables)
        });

    //println!("done");

    system.varmap = new_varmap;
    reductions
}

fn apply_all_substitutions<N: Numeric>(
    mut equation: EquationViewMut<N>,
    substitutions: &[Option<Equation<N>>],
    storage_variables: usize,
) {
    let mut new_equation = Equation {
        data: vec![N::from(0); storage_variables + 1],
    };
    let mut scratch = N::from(0);

    let mut new_equation = new_equation.as_mut();

    new_equation.get_result().clone_from(equation.get_result());

    for (i, s) in substitutions.iter().enumerate() {
        if s.is_none() {
            new_equation
                .get_coefficient(i)
                .clone_from(equation.get_coefficient(i));
        } else {
            new_equation.get_coefficient(i).assign(0);
        }
    }

    //println!("Applying all substitutions:");

    for (subst_idx, subst) in substitutions.iter().enumerate() {
        if let Some(subst) = subst {
            //println!("current result {}", new_equation.display(varmap));
            let factor = equation.get_coefficient(subst_idx);

            for (data_idx, subst_coef) in subst.data.iter().enumerate() {
                let s = &mut scratch;
                s.clone_from(factor);
                *s *= subst_coef;

                new_equation.data[data_idx] += &*s
            }
        }
    }

    //println!("final result {}", new_equation.display(varmap));

    equation.copy_from(new_equation.as_ref());
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
