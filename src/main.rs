use std::path::PathBuf;

use structopt::StructOpt;

mod parser;

pub struct EquationSystem {
    variables: usize,
    equations: usize,
    data: Vec<i32>,
}

#[derive(Debug, StructOpt)]
struct Config {
    input: PathBuf
}

fn main() -> anyhow::Result<()> {
    let config = Config::from_args();

    let input = std::fs::read_to_string(config.input)?;

    let result = parser::parse(&input)?;

    Ok(())
}
