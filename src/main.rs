use std::path::PathBuf;

use structopt::StructOpt;

mod parser;
mod system;
mod algo;
mod math;
mod util;

#[derive(Debug, StructOpt)]
struct Config {
    input: PathBuf
}

fn main() -> anyhow::Result<()> {
    let config = Config::from_args();

    let input = std::fs::read_to_string(config.input)?;

    let mut result = parser::parse(&input)?;

    println!("{}", result.equations_display());

    let result = algo::solve_equation(&mut result);

    if result == algo::Result::Sat {
        println!("Solvable");
    } else {
        println!("Unsolvable");
    }

    Ok(())
}
