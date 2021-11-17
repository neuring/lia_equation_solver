use std::{convert::TryInto, ops::Not};

use thiserror::Error;

use crate::system::System;

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("Value `{0}` is not a number.")]
    NotANumber(String),

    #[error("Missing value {0}.")]
    MissingValue(&'static str),

    #[error("Unexpected end of input.")]
    UnexpectedEndOfInput,

    #[error("Index {0} too large for specified number of variables.")]
    IndexTooLarge(usize),

    #[error("Expected trailing zero is not zero ({0})..")]
    TrailingZeroNotZero(i32),

    #[error("Index {0} is invalid")]
    InvalidIndex(i32)
}

use ParseError::*;

pub type Result<T> = std::result::Result<T, ParseError>;

struct HeaderData {
    equations: usize,
    variables: usize,
}

fn parse_header(header: &str) -> Result<HeaderData> {
    let mut values = header.split_whitespace();

    let number_equations =
        values.next().ok_or(MissingValue("number of equations."))?;
    let number_equations = number_equations
        .parse::<usize>()
        .map_err(|_| NotANumber(number_equations.to_string()))?;

    let number_variables =
        values.next().ok_or(MissingValue("number of variables."))?;
    let number_variables = number_variables
        .parse::<usize>()
        .map_err(|_| NotANumber(number_variables.to_string()))?;

    Ok(HeaderData {
        equations: number_equations,
        variables: number_variables,
    })
}

fn parse_equation<'a>(input: &'a str, storage: &mut [i32]) -> Result<()> {
    let equation_size = storage.len();

    let mut values = input
        .split_whitespace()
        .map(|val| val.parse::<i32>().map_err(|_| NotANumber(val.to_string())));

    let number_terms = values.next().ok_or(MissingValue("number of terms"))??;

    for _ in 0..number_terms - 1 {
        let coefficient = values.next().ok_or(MissingValue("coefficient"))??;
        let index = values.next().ok_or(MissingValue("variable index"))??;

        let index: usize = index.try_into().map_err(|_| InvalidIndex(index))?;

        if index < equation_size {
            storage[index] = coefficient;
        } else {
            return Err(IndexTooLarge(index));
        }
    }

    let equation_result = values.next().ok_or(MissingValue("equation result"))??;

    storage[0] = equation_result;

    let expected_zero = values.next().ok_or(MissingValue("trailing zero"))??;

    if expected_zero != 0 {
        return Err(TrailingZeroNotZero(expected_zero));
    }

    Ok(())
}

fn is_empty_line(line: &str) -> Option<&str> {
    let line = line.trim();

    (line.is_empty() || line.starts_with("#")).not().then(|| line)
}

pub fn parse(input: &str) -> Result<System> {
    let mut lines = input.lines().filter_map(is_empty_line);

    let header = parse_header(lines.next().ok_or(UnexpectedEndOfInput)?)?;

    let equation_size = header.variables + 1;

    let mut result = vec![0; header.equations * equation_size];

    for equation_idx in 0..header.equations {
        let buffer_start = equation_idx * equation_size;
        let equation_buffer =
            &mut result[buffer_start..buffer_start + equation_size];

        parse_equation(
            lines.next().ok_or(UnexpectedEndOfInput)?,
            equation_buffer,
        )?;
    }

    Ok(System::new(
        header.variables,
        header.equations,
        result,
    ))
}
