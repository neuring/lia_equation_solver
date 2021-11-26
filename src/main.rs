use std::path::PathBuf;

use structopt::StructOpt;

use crate::system::{VariableIndex, System};

mod algo;
mod math;
mod parser;
mod system;
mod util;
mod numeric;

#[derive(Debug, StructOpt)]
struct Config {
    #[structopt(short, long)]
    dump_dot: Option<PathBuf>,

    #[structopt(short, long)]
    verify: bool,

    input: PathBuf,
}

type N = i64;

fn main() -> anyhow::Result<()> {
    let config = Config::from_args();

    let input = std::fs::read_to_string(config.input)?;

    let mut system: System<N> = parser::parse(&input)?;

    let original_system = system.clone();

    println!("{}", system.equations_display());

    let result = algo::solve_equation(&mut system);

    if result == algo::Result::Sat {
        println!("Solvable");

        let result = system
            .reconstruction
            .evaluate_with_zeroes(system.next_var_index, &mut N::from(0));

        for (i, res) in result.iter().copied().enumerate() {
            if i >= system.storage.variables {
                break;
            }
            if let Some(res) = res {
                println!(
                    "{} = {}",
                    util::fmt_variable(VariableIndex(i), system.storage.variables),
                    res
                )
            }
        }

        if config.verify {
            let assignment: Vec<_> =
                result.iter().map(|x| x.unwrap_or(-123)).collect();

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

    if let Some(dump_path) = config.dump_dot {
        let f = std::fs::File::create(dump_path)?;
        system
            .reconstruction
            .dump_dot(system.storage.variables, f)?;
    }

    Ok(())
}
