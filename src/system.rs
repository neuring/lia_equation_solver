use std::{collections::HashMap, fmt};

use itertools;

use crate::util;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VariableIndex(pub usize);

#[derive(Debug, Clone)]
pub struct System {
    pub varmap: Vec<VariableIndex>,
    pub next_var_index: usize,

    pub reconstruction: Reconstruction,

    pub scratch_pad: Vec<i64>,

    pub storage: EquationStorage,
}

impl System {
    pub fn new(variables: usize) -> Self {
        Self {
            varmap: (0..variables).map(|i| VariableIndex(i)).collect(),
            next_var_index: variables,
            scratch_pad: vec![0; variables],

            reconstruction: Reconstruction::new(),

            storage: EquationStorage::new(variables),
        }
    }

    pub fn new_variable(&mut self) -> VariableIndex {
        let var = VariableIndex(self.next_var_index);
        self.next_var_index += 1;
        var
    }

    pub fn map_variable(&mut self, term_idx: usize, variable_idx: VariableIndex) {
        self.varmap[term_idx] = variable_idx;
    }

    pub fn add_equation(&mut self) -> EquationViewMut {
        self.storage.add_equation()
    }

    pub fn get_storage(&self) -> &EquationStorage {
        &self.storage
    }

    pub fn get_storage_mut(&mut self) -> &mut EquationStorage {
        &mut self.storage
    }

    pub fn equations_display(&self) -> impl fmt::Display + '_ {
        struct DisplayableEquations<'a>(&'a System);

        impl<'a> fmt::Display for DisplayableEquations<'a> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                for equation in
                    self.0.storage.iter_equations().filter(|eq| !eq.is_empty())
                {
                    writeln!(f, "{}", equation.equation_display(&self.0.varmap))?;
                }
                Ok(())
            }
        }

        DisplayableEquations(self)
    }

    /// Evaluates if the assignment is a solution to the system of equations.
    /// `assignment` maps variable index to its assignment.
    pub fn evaluate(&self, assignment: &[i64]) -> Result<(), EquationView> {
        for equation in self.storage.iter_equations() {
            let evaluate_result: i64 = self
                .varmap
                .iter()
                .copied()
                .zip(equation.iter_coefficients())
                .map(|(variable_idx, coeff)| assignment[variable_idx.0] * coeff)
                .sum();

            if evaluate_result != equation.get_result() {
                return Err(equation);
            }
        }

        Ok(())
    }
}

#[derive(Clone)]
pub struct EquationStorage {
    pub variables: usize,
    pub equations: usize,

    pub data: Vec<i64>,
}

impl fmt::Debug for EquationStorage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let equations: Vec<_> = self.iter_equations().collect();

        f.debug_struct("System")
            .field("variables", &self.variables)
            .field("equations", &self.equations)
            .field("data", &equations)
            .finish()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct EquationView<'a> {
    data: &'a [i64],
}

#[derive(Debug)]
pub struct EquationViewMut<'a> {
    pub data: &'a mut [i64],
}

#[derive(Debug, Clone)]
pub struct Equation {
    data: Vec<i64>,
}

impl EquationStorage {
    pub fn new(variables: usize) -> Self {
        Self {
            variables,
            equations: 0,
            data: Vec::new(),
        }
    }

    pub fn add_equation(&mut self) -> EquationViewMut {
        self.data
            .resize(self.data.len() + self.get_equation_size(), 0);

        let equation_idx = self.equations;

        self.equations += 1;

        let view = self.get_equation_mut(equation_idx);

        view
    }

    /// Returns the number of terms in an equation.
    pub fn get_terms(&self) -> usize {
        return self.variables;
    }

    /// Returns the number of equations in the system of equations.
    pub fn get_equations(&self) -> usize {
        return self.variables;
    }

    pub fn get_equation(&self, idx: usize) -> EquationView {
        assert!(idx < self.get_equations());

        let start = self.get_equation_size() * idx;
        let data = &self.data[start..start + self.get_equation_size()];

        EquationView { data }
    }

    pub fn get_equation_mut(&mut self, idx: usize) -> EquationViewMut {
        assert!(idx < self.get_equations());

        let size = self.get_equation_size();
        let start = self.get_equation_size() * idx;
        let data = &mut self.data[start..start + size];

        EquationViewMut { data }
    }

    /// Returns the number of integers an equation contains.
    /// It is the number of terms plus one, to account for the equation result.
    fn get_equation_size(&self) -> usize {
        return self.variables + 1;
    }

    pub fn iter_equations(&self) -> impl Iterator<Item = EquationView<'_>> + '_ {
        self.data
            .chunks_exact(self.get_equation_size())
            .map(|equation_data| EquationView {
                data: equation_data,
            })
    }

    pub fn iter_equations_mut(
        &mut self,
    ) -> impl Iterator<Item = EquationViewMut<'_>> + '_ {
        let equation_size = self.get_equation_size();
        self.data
            .chunks_exact_mut(equation_size)
            .map(|equation_data| EquationViewMut {
                data: equation_data,
            })
    }
}

impl<'a> EquationView<'a> {
    pub fn get_result(&self) -> i64 {
        return self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&self, idx: usize) -> i64 {
        return self.data[idx + 1];
    }

    pub fn iter_coefficients(&self) -> impl Iterator<Item = i64> + '_ {
        self.data.iter().skip(1).copied()
    }

    pub fn to_owned(&self) -> Equation {
        Equation {
            data: self.data.to_owned(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.data.iter().all(|&x| x == 0)
    }

    pub fn equation_display(
        &self,
        varmap: &'a Vec<VariableIndex>,
    ) -> impl fmt::Display + 'a {
        struct EquationDisplay<'a> {
            equation: EquationView<'a>,
            varmap: &'a Vec<VariableIndex>,
        }

        impl<'a> fmt::Display for EquationDisplay<'a> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                let mut terms: Vec<_> = self
                    .equation
                    .iter_coefficients()
                    .enumerate()
                    .filter(|&(_, coef)| coef != 0)
                    .map(|(idx, coef)| (coef, self.varmap[idx]))
                    .collect();

                terms.sort_by_key(|(_, i)| i.0);

                let equation_lhs = if terms.is_empty() {
                    "0".to_owned()
                } else {
                    itertools::join(
                        terms.iter().map(|&(coef, idx)| {
                            format!(
                                "{}{}",
                                coef,
                                util::fmt_variable(idx, self.varmap.len()),
                            )
                        }),
                        " + ",
                    )
                };

                write!(f, "{} = {}", equation_lhs, self.equation.get_result())
            }
        }

        EquationDisplay {
            equation: *self,
            varmap,
        }
    }
}

impl<'a> EquationViewMut<'a> {
    pub fn get_result(&mut self) -> &mut i64 {
        return &mut self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&mut self, idx: usize) -> &mut i64 {
        return &mut self.data[idx + 1];
    }

    pub fn iter_coefficients(&mut self) -> impl Iterator<Item = &mut i64> + '_ {
        self.data.iter_mut().skip(1)
    }

    pub fn get_coefficient_slice(&mut self) -> &mut [i64] {
        &mut self.data[1..]
    }

    /// Sets all coefficient to zero.
    pub fn clear(&mut self) {
        self.data.fill(0);
    }

    pub fn to_owned(&self) -> Equation {
        Equation {
            data: self.data.to_owned(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.data.iter().all(|&x| x == 0)
    }

    pub fn display(
        &'a self,
        varmap: &'a Vec<VariableIndex>,
    ) -> impl fmt::Display + 'a {
        EquationView { data: &self.data }.equation_display(varmap)
    }
}

impl Equation {
    pub fn get_result_mut(&mut self) -> &mut i64 {
        return &mut self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient_mut(&mut self, idx: usize) -> &mut i64 {
        return &mut self.data[idx + 1];
    }

    pub fn iter_coefficients_mut(&mut self) -> impl Iterator<Item = &mut i64> + '_ {
        self.data.iter_mut().skip(1)
    }

    pub fn get_result(&self) -> i64 {
        return self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&self, idx: usize) -> i64 {
        return self.data[idx + 1];
    }

    pub fn iter_coefficients(&self) -> impl Iterator<Item = i64> + '_ {
        self.data.iter().skip(1).copied()
    }

    pub fn display<'a>(
        &'a self,
        varmap: &'a Vec<VariableIndex>,
    ) -> impl fmt::Display + '_ {
        EquationView { data: &self.data }.equation_display(varmap)
    }
}

#[derive(Debug, Clone)]
pub struct Reconstruction {
    tree: HashMap<VariableIndex, ReconstructedEquation>,
}

impl Reconstruction {
    pub fn new() -> Self {
        Self {
            tree: HashMap::new(),
        }
    }

    pub fn add(
        &mut self,
        var: VariableIndex,
        terms: Vec<(VariableIndex, i64)>,
        constant: i64,
    ) {
        assert!(!self.tree.contains_key(&var));

        self.tree
            .insert(var, ReconstructedEquation { constant, terms });
    }

    pub fn dump_dot(
        &self,
        total_variables: usize,
        mut w: impl std::io::Write,
    ) -> Result<(), impl std::error::Error> {
        writeln!(w, "digraph {{")?;

        for (var, def) in self.tree.iter() {
            let fmt_var = util::fmt_variable(*var, total_variables);
            writeln!(
                w,
                "{} [label=\"{}: {}\" ordering=\"out\"];",
                fmt_var, fmt_var, def.constant,
            )?;
            for (child_var, coef) in &def.terms {
                writeln!(
                    w,
                    "{} -> {} [label={}];",
                    fmt_var,
                    util::fmt_variable(*child_var, total_variables),
                    coef
                )?;
            }
        }

        writeln!(w, "}}")
    }

    pub fn evaluate_with_zeroes(&self, total_variables: usize) -> Vec<Option<i64>> {
        let mut visited = vec![None; total_variables];

        for var in self.tree.keys() {
            self.evaluate_recursive(*var, &mut visited);
        }

        visited
    }

    fn evaluate_recursive(
        &self,
        current: VariableIndex,
        visited: &mut Vec<Option<i64>>,
    ) {
        if visited[current.0].is_some() {
            return;
        }

        let def = self.tree.get(&current);

        if let Some(def) = def {
            let mut value = def.constant;

            for &(child, coef) in &def.terms {
                if visited[child.0].is_none() {
                    self.evaluate_recursive(child, visited);
                }

                value += visited[child.0].unwrap() * coef;
            }

            visited[current.0] = Some(value);
        } else {
            visited[current.0] = Some(0);
        }
    }
}

#[derive(Debug, Clone)]
struct ReconstructedEquation {
    constant: i64,
    terms: Vec<(VariableIndex, i64)>,
}
