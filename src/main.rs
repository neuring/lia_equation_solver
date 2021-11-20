use std::path::PathBuf;

use structopt::StructOpt;

use crate::system::VariableIndex;

mod algo;
mod math;
mod parser;
mod system;
mod util;

#[derive(Debug, StructOpt)]
struct Config {
    #[structopt(short, long)]
    dump_dot: Option<PathBuf>,

    input: PathBuf,
}

fn main() -> anyhow::Result<()> {
    let config = Config::from_args();

    let input = std::fs::read_to_string(config.input)?;

    let mut system = parser::parse(&input)?;

    println!("{}", system.equations_display());

    let result = algo::solve_equation(&mut system);

    if result == algo::Result::Sat {
        println!("Solvable");

        let result = system
            .reconstruction
            .evaluate_with_zeroes(system.next_var_index);

        for (i, res) in result.into_iter().enumerate() {
            if i >= system.storage.variables { break; }
            if let Some(res) = res {
                println!(
                    "{} = {}",
                    util::fmt_variable(VariableIndex(i), system.storage.variables),
                    res
                )
            }
        }
    } else {
        println!("Unsolvable");
    }

    if let Some(dump_path) = config.dump_dot {
        let f = std::fs::File::create(dump_path)?;
        system.reconstruction.dump_dot(system.storage.variables, f)?;
    }

    Ok(())
}
