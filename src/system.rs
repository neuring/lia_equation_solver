use std::{collections::HashMap, fmt};

use itertools;

use crate::{numeric::Numeric, util};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VariableIndex(pub usize);

#[derive(Debug, Clone)]
pub struct System<N> {
    pub varmap: Vec<VariableIndex>,

    pub alive_terms: Vec<bool>,

    pub killed_variables: usize,

    pub next_var_index: usize,

    pub starting_variables: usize,

    pub reconstruction: Reconstruction<N>,

    pub storage: EquationStorage<N>,
}

impl<N: Numeric> System<N> {
    pub fn new(variables: usize) -> Self {
        Self {
            varmap: (0..variables).map(|i| VariableIndex(i)).collect(),
            next_var_index: variables,

            alive_terms: vec![true; variables],
            killed_variables: 0,

            starting_variables: variables,

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

    pub fn kill_variable(&mut self, term_idx: usize) {
        self.killed_variables += 1;
        self.alive_terms[term_idx] = false;
    }

    pub fn add_equation(&mut self) -> EquationViewMut<'_, N> {
        self.storage.add_equation()
    }

    pub fn get_storage(&self) -> &EquationStorage<N> {
        &self.storage
    }

    pub fn equations_display(&self) -> impl fmt::Display + '_ {
        struct DisplayableEquations<'a, N>(&'a System<N>);

        impl<'a, N: Numeric> fmt::Display for DisplayableEquations<'a, N> {
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
    pub fn evaluate(&self, assignment: &[N]) -> Result<(), EquationView<'_, N>> {
        for equation in self.storage.iter_equations() {
            //.map(|(variable_idx, coeff)| assignment[variable_idx.0] * coeff)
            let mut scratch = N::from(0);
            let mut sum = N::from(0);
            self.varmap
                .iter()
                .copied()
                .zip(equation.iter_coefficients())
                .for_each(|(variable_idx, coeff)| {
                    scratch.clone_from(&assignment[variable_idx.0]);
                    scratch *= coeff;
                    sum += &scratch;
                });

            sum.negate();

            if &sum != equation.get_result() {
                return Err(equation);
            }
        }

        Ok(())
    }

    pub fn resize(&mut self) {
        let new_variables = self.starting_variables - self.killed_variables;
        let new_equation_size = new_variables + 1;

        let old_equation_size = self.storage.get_equation_size();

        let data = &mut self.storage.data;

        let equations = self.storage.equations;

        let mut empty_equations = 0;

        let mut dst_eq_idx = 0;
        for src_eq_idx in 0..equations {
            let src_eq_start = src_eq_idx * old_equation_size;
            let eq_data = &data[src_eq_start..src_eq_start + old_equation_size];

            // Is src equation empty? -> continue
            if eq_data.iter().all(|i| i.cmp_zero().is_eq()) {
                empty_equations += 1;
                continue;
            }

            // swap elements of src equation with smaller dst equation
            let dst_eq_start = dst_eq_idx * new_equation_size;

            // swap equation result values
            data.swap(dst_eq_start, src_eq_start);

            // swap equation coefficients
            let mut dst_i = 1;
            for src_i in 1..old_equation_size {
                // in data, coefficients start at one, in varmap and alive_terms at zero
                if self.alive_terms[src_i - 1] {
                    data.swap(dst_eq_start + dst_i, src_eq_start + src_i);
                    dst_i += 1;
                }
            }
            assert_eq!(dst_i, new_equation_size);

            dst_eq_idx += 1;
        }

        // resize varmap and alive_terms for new smaller number of active variables
        assert_eq!(self.varmap.len(), self.alive_terms.len());

        let mut alive = self.alive_terms.iter().copied();
        self.varmap.retain(|_| alive.next().unwrap());
        self.alive_terms.resize(self.varmap.len(), true);
        self.alive_terms.fill(true);

        assert_eq!(self.alive_terms.len(), new_variables);

        self.storage.equations -= empty_equations;
        self.storage.variables = new_variables;
        data.truncate(self.storage.equations * (new_variables + 1));
        data.shrink_to_fit()
    }
}

#[derive(Clone)]
pub struct EquationStorage<N> {
    pub variables: usize,
    pub equations: usize,

    pub data: Vec<N>,
}

impl<N: fmt::Debug> fmt::Debug for EquationStorage<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let equations: Vec<_> = self.iter_equations().collect();

        f.debug_struct("System")
            .field("variables", &self.variables)
            .field("equations", &self.equations)
            .field("data", &equations)
            .finish()
    }
}

#[derive(Debug)]
pub struct EquationView<'a, N> {
    data: &'a [N],
}

impl<'a, N> Clone for EquationView<'a, N> {
    fn clone(&self) -> Self {
        Self { data: self.data }
    }
}

impl<'a, N> Copy for EquationView<'a, N> {}

#[derive(Debug)]
pub struct EquationViewMut<'a, N> {
    pub data: &'a mut [N],
}

#[derive(Debug, Clone)]
pub struct Equation<N> {
    data: Vec<N>,
}

impl<N> EquationStorage<N> {
    pub fn new(variables: usize) -> Self {
        Self {
            variables,
            equations: 0,
            data: Vec::new(),
        }
    }

    /// Returns the number of equations in the system of equations.
    pub fn get_equations(&self) -> usize {
        return self.variables;
    }

    pub fn get_equation(&self, idx: usize) -> EquationView<'_, N> {
        assert!(idx < self.get_equations());

        let start = self.get_equation_size() * idx;
        let data = &self.data[start..start + self.get_equation_size()];

        EquationView { data }
    }

    pub fn get_equation_mut(&mut self, idx: usize) -> EquationViewMut<'_, N> {
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

    pub fn iter_equations(&self) -> impl Iterator<Item = EquationView<'_, N>> + '_ {
        self.data
            .chunks_exact(self.get_equation_size())
            .take(self.equations)
            .map(|equation_data| EquationView {
                data: equation_data,
            })
    }

    pub fn iter_equations_mut(
        &mut self,
    ) -> impl Iterator<Item = EquationViewMut<'_, N>> + '_ {
        let equation_size = self.get_equation_size();
        self.data
            .chunks_exact_mut(equation_size)
            .map(|equation_data| EquationViewMut {
                data: equation_data,
            })
    }
}

impl<N: Numeric> EquationStorage<N> {
    pub fn add_equation(&mut self) -> EquationViewMut<'_, N> {
        self.data
            .resize(self.data.len() + self.get_equation_size(), N::from(0));

        let equation_idx = self.equations;

        self.equations += 1;

        let view = self.get_equation_mut(equation_idx);

        view
    }

    pub fn print_value_stats(&self, idx: usize, mut writer: impl std::io::Write) {
        self.data
            .iter()
            .filter(|i| i.cmp_zero().is_ne())
            .for_each(|i| writeln!(writer, "{} {}", idx, i).unwrap());
    }
}

impl<'a, N> EquationView<'a, N> {
    pub fn get_result(&self) -> &N {
        return &self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    #[allow(unused)]
    pub fn get_coefficient(&self, idx: usize) -> &'a N {
        return &self.data[idx + 1];
    }

    pub fn iter_coefficients(&self) -> impl Iterator<Item = &N> + '_ {
        self.data.iter().skip(1)
    }
}

#[allow(unused)]
impl<'a, N: Clone> EquationView<'a, N> {
    pub fn to_owned(&self) -> Equation<N> {
        Equation {
            data: self.data.to_owned(),
        }
    }
}

impl<'a, N: Numeric> EquationView<'a, N> {
    pub fn is_empty(&self) -> bool {
        self.data.iter().all(|x| x.cmp_zero().is_eq())
    }

    pub fn equation_display(
        &self,
        varmap: &'a Vec<VariableIndex>,
    ) -> impl fmt::Display + 'a {
        struct EquationDisplay<'a, N> {
            equation: EquationView<'a, N>,
            varmap: &'a Vec<VariableIndex>,
        }

        impl<'a, N: Numeric> fmt::Display for EquationDisplay<'a, N> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                let mut terms: Vec<_> = self
                    .equation
                    .iter_coefficients()
                    .enumerate()
                    .filter(|&(_, coef)| coef.cmp_zero().is_ne())
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

                write!(f, "{} + {}", self.equation.get_result(), equation_lhs)
            }
        }

        EquationDisplay {
            equation: *self,
            varmap,
        }
    }
}

impl<'a, N: Numeric> EquationViewMut<'a, N> {
    pub fn get_result(&mut self) -> &mut N {
        return &mut self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&mut self, idx: usize) -> &mut N {
        return &mut self.data[idx + 1];
    }

    pub fn iter_coefficients(&mut self) -> impl Iterator<Item = &mut N> + '_ {
        self.data.iter_mut().skip(1)
    }

    pub fn get_coefficient_slice(&mut self) -> &mut [N] {
        &mut self.data[1..]
    }

    /// Sets all coefficient to zero.
    pub fn clear(&mut self) {
        self.data.iter_mut().for_each(|i| i.assign(0));
    }

    pub fn to_owned(&self) -> Equation<N> {
        Equation {
            data: self.data.to_owned(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.data.iter().all(|x| x.cmp_zero().is_eq())
    }

    #[allow(unused)]
    pub fn display(
        &'a self,
        varmap: &'a Vec<VariableIndex>,
    ) -> impl fmt::Display + 'a {
        EquationView { data: &self.data }.equation_display(varmap)
    }

    pub fn copy_into(&mut self, other: EquationView<'_, N>) {
        for (s, o) in self.data.iter_mut().zip(other.data.iter()) {
            s.clone_from(o);
        }
    }

    pub fn into_ref(self) -> EquationView<'a, N> {
        EquationView { data: &*self.data }
    }
}

impl<N: Numeric> Equation<N> {
    pub fn get_result(&self) -> &N {
        return &self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&self, idx: usize) -> &N {
        return &self.data[idx + 1];
    }

    pub fn iter_coefficients(&self) -> impl Iterator<Item = &N> + '_ {
        self.data.iter().skip(1)
    }

    #[allow(unused)]
    pub fn display<'a>(
        &'a self,
        varmap: &'a Vec<VariableIndex>,
    ) -> impl fmt::Display + '_ {
        EquationView { data: &self.data }.equation_display(varmap)
    }
}

#[derive(Debug, Clone)]
pub struct Reconstruction<N> {
    tree: HashMap<VariableIndex, ReconstructedEquation<N>>,
}

impl<N> Reconstruction<N> {
    pub fn new() -> Self {
        Self {
            tree: HashMap::new(),
        }
    }

    pub fn add(
        &mut self,
        var: VariableIndex,
        terms: Vec<(VariableIndex, N)>,
        constant: N,
    ) {
        assert!(!self.tree.contains_key(&var));

        self.tree
            .insert(var, ReconstructedEquation { constant, terms });
    }
}

impl<N: fmt::Display> Reconstruction<N> {
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
}

impl<N: Numeric> Reconstruction<N> {
    pub fn evaluate_with_zeroes(
        &self,
        total_variables: usize,
        scratch: &mut N,
    ) -> Vec<Option<N>> {
        let mut visited = vec![None; total_variables];

        for var in self.tree.keys() {
            self.evaluate_recursive(*var, &mut visited, scratch);
        }

        visited
    }

    fn evaluate_recursive(
        &self,
        current: VariableIndex,
        visited: &mut Vec<Option<N>>,
        scratch: &mut N,
    ) {
        if visited[current.0].is_some() {
            return;
        }

        let def = self.tree.get(&current);

        if let Some(def) = def {
            let mut value = def.constant.clone();

            for (child, coef) in &def.terms {
                if visited[child.0].is_none() {
                    self.evaluate_recursive(*child, visited, scratch);
                }

                scratch.clone_from(visited[child.0].as_ref().unwrap());
                *scratch *= coef;
                value += &*scratch;
            }

            visited[current.0] = Some(value);
        } else {
            visited[current.0] = Some(N::from(0));
        }
    }
}

#[derive(Debug, Clone)]
struct ReconstructedEquation<N> {
    constant: N,
    terms: Vec<(VariableIndex, N)>,
}
