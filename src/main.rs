use std::path::PathBuf;

use anyhow::Context;
use structopt::StructOpt;

use crate::system::{System, VariableIndex};

mod algo;
mod math;
mod numeric;
mod parser;
mod system;
mod util;

#[derive(Debug, StructOpt)]
struct Config {
    #[structopt(subcommand)]
    command: Command,

    /// Path to input file
    input: PathBuf,
}

#[derive(Debug, StructOpt)]
enum Command {
    /// Convert input to smt-lib format
    ToSMTLib { output: Option<PathBuf> },

    /// Solve system of equation.
    Solve(SolverConfig),
}

#[derive(Debug, StructOpt)]
struct SolverConfig {
    /// dump solution as dot graph
    #[structopt(short, long)]
    dump_dot: Option<PathBuf>,

    /// Do not verify the solution
    #[structopt(short, long)]
    no_verify: bool,
}

type N = rug::Integer;

fn solver_main(mut system: System<N>, config: &SolverConfig) -> anyhow::Result<()> {
    let original_system = system.clone();

    if system.starting_variables <= 10 {
        println!("{}", system.equations_display());
    }

    let result = algo::solve_equation(&mut system);

    if result == algo::Result::Sat {
        println!("Solvable");

        if let Some(dump_path) = config.dump_dot.as_ref() {
            let f = std::fs::File::create(dump_path)?;
            system
                .reconstruction
                .dump_dot(system.storage.variables, f)?;
        }

        let result = system
            .reconstruction
            .evaluate_solution(system.var_generator.next_var_index, &mut N::from(0));

        for (i, res) in result.iter().cloned().enumerate() {
            if i >= system.starting_variables {
                break;
            }
            if let Some(res) = res {
                println!(
                    "{} = {}",
                    util::fmt_variable(VariableIndex(i), system.starting_variables),
                    res.display_solution()
                )
            }
        }

        if !config.no_verify {
            let assignment: Vec<_> = result
                .into_iter()
                .map(|x| x.map(|e| e.constant).unwrap_or_else(|| N::from(0)))
                .collect();

            match original_system.evaluate(&assignment) {
                Ok(()) => println!("solution verified."),
                Err(equation) => println!(
                    "wrong solution: {}",
                    equation.equation_display(&original_system.varmap)
                ),
            }
        }
    } else {
        println!("Unsolvable");
    }

    if result == algo::Result::Sat {
        Ok(())
    } else {
        Err(anyhow::anyhow!(""))
    }
}

fn smtlib_main(system: System<N>, out_path: Option<PathBuf>) -> anyhow::Result<()> {
    todo!()
}

fn main() -> anyhow::Result<()> {
    let config = Config::from_args();

    let input = std::fs::read_to_string(&config.input).with_context(|| {
        format!("Couldn't read file '{}'", &config.input.to_string_lossy())
    })?;

    let system: System<N> =
        parser::parse(&input).context("Failed to parse input.")?;

    match config.command {
        Command::ToSMTLib { output } => todo!(),
        Command::Solve(solver_config) => solver_main(system, &solver_config),
    }
}
