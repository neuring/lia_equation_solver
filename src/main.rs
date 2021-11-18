use std::path::PathBuf;

use structopt::StructOpt;

mod parser;
mod system;
mod algo;
mod math;

#[derive(Debug, StructOpt)]
struct Config {
    input: PathBuf
}

fn main() -> anyhow::Result<()> {
    let config = Config::from_args();

    let input = std::fs::read_to_string(config.input)?;

    let mut result = parser::parse(&input)?;

    dbg!(&result);

    algo::solve_equation(&mut result);

    dbg!(&result);

    Ok(())
}
