use std::fmt;

pub struct System {
    variables: usize,
    equations: usize,
    data: Vec<i32>,
}

impl fmt::Debug for System {
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
pub struct EquationView<'a> {
    data: &'a [i32],
}

#[derive(Debug)]
pub struct EquationViewMut<'a> {
    data: &'a mut [i32],
}

#[derive(Debug)]
pub struct Equation {
    data: Vec<i32>,
}

impl System {
    pub fn new(variables: usize, equations: usize, data: Vec<i32>) -> Self {
        assert_eq!((variables + 1) * equations, data.len());

        Self {
            variables,
            equations,
            data,
        }
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
    pub fn get_result(&self) -> i32 {
        return self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&self, idx: usize) -> i32 {
        return self.data[idx + 1];
    }

    pub fn iter_coefficients(&self) -> impl Iterator<Item = i32> + '_ {
        self.data.iter().skip(1).copied()
    }

    pub fn to_owned(&self) -> Equation {
        Equation {
            data: self.data.to_owned(),
        }
    }
}

impl<'a> EquationViewMut<'a> {
    pub fn get_result(&mut self) -> &mut i32 {
        return &mut self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&mut self, idx: usize) -> &mut i32 {
        return &mut self.data[idx + 1];
    }

    pub fn iter_coefficients(&mut self) -> impl Iterator<Item = &mut i32> + '_ {
        self.data.iter_mut().skip(1)
    }

    pub fn get_coefficient_slice(&mut self) -> &mut [i32] {
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
}

impl Equation {
    pub fn get_result_mut(&mut self) -> &mut i32 {
        return &mut self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient_mut(&mut self, idx: usize) -> &mut i32 {
        return &mut self.data[idx + 1];
    }

    pub fn iter_coefficients_mut(&mut self) -> impl Iterator<Item = &mut i32> + '_ {
        self.data.iter_mut().skip(1)
    }

    pub fn get_result(&self) -> i32 {
        return self.data[0];
    }

    /// Returns a coefficient of the equation.
    /// `idx` is not the variable index, but the index where it is stored in memory.
    pub fn get_coefficient(&self, idx: usize) -> i32 {
        return self.data[idx + 1];
    }

    pub fn iter_coefficients(&self) -> impl Iterator<Item = i32> + '_ {
        self.data.iter().skip(1).copied()
    }
}
