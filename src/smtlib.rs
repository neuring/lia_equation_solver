use std::io::Write;

use crate::{
    numeric::Numeric,
    system::{EquationView, System},
};

fn write_prolog<N: Numeric>(
    system: &System<N>,
    out: &mut dyn Write,
) -> anyhow::Result<()> {
    writeln!(out, "(set-logic QF_LIA)")?;

    for i in 0..system.starting_variables {
        writeln!(out, "(declare-const x{} Int)", i)?;
    }

    Ok(())
}

fn write_epilog(out: &mut dyn Write) -> anyhow::Result<()> {
    writeln!(out, "(check-sat)")?;
    writeln!(out, "(exit)")?;
    Ok(())
}

fn write_equation<N: Numeric>(
    equation: EquationView<N>,
    out: &mut dyn Write,
) -> anyhow::Result<()> {
    let mut result = N::from(-1);
    result *= equation.get_result();

    let terms = equation
        .iter_coefficients()
        .enumerate()
        .filter(|(_, c)| c.cmp_zero().is_ne())
        .map(|(idx, c)| format!("(* {} x{})", c, idx));

    let terms = itertools::join(terms, " ");

    writeln!(out, "(assert (= {} (+ {})))", result, terms)?;

    Ok(())
}

fn write_body<N: Numeric>(
    system: &System<N>,
    out: &mut dyn Write,
) -> anyhow::Result<()> {
    for equation in system.storage.iter_equations() {
        write_equation(equation, out)?;
    }
    Ok(())
}

pub fn write_smtlib<N: Numeric>(
    system: &System<N>,
    out: &mut dyn Write,
) -> anyhow::Result<()> {
    write_prolog(system, out)?;

    write_body(&system, out)?;

    write_epilog(out)?;

    Ok(())
}
