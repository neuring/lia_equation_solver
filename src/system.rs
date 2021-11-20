use std::fmt;

use itertools;

use crate::util;

#[derive(Debug, Clone, Copy)]
pub struct VariableIndex(usize);

pub struct System {
    varmap: Vec<VariableIndex>,
    next_var_index: usize,



    pub scratch_pad: Vec<i64>,

    pub storage: EquationStorage,
}

impl System {
    pub fn new(variables: usize) -> Self {
        Self {
            varmap: (0..variables).map(|i| VariableIndex(i)).collect(),
            next_var_index: variables,
            scratch_pad: vec![0; variables],

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
                    let mut terms: Vec<_> = equation
                        .iter_coefficients()
                        .enumerate()
                        .filter(|&(_, coef)| coef != 0)
                        .map(|(idx, coef)| (coef, self.0.varmap[idx]))
                        .collect();

                    terms.sort_by_key(|(_, i)| i.0);

                    let equation_lhs = if terms.is_empty() {
                        "0".to_owned()
                    } else {
                        itertools::join(
                            terms.iter().map(|&(coef, idx)| {
                                let var = if idx.0 < self.0.storage.variables {
                                    'x'
                                } else {
                                    'y'
                                };

                                format!(
                                    "{}{}{}",
                                    coef,
                                    var,
                                    util::subscript(idx.0 as u32)
                                )
                            }),
                            " + ",
                        )
                    };

                    writeln!(f, "{} = {}", equation_lhs, equation.get_result())?;
                }
                Ok(())
            }
        }

        DisplayableEquations(self)
    }
}

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
    data: &'a mut [i64],
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
}
