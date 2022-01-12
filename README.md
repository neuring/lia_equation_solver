# Linear Integer Arithmetic Equation solver

## Build
To build this solver it is recommended to use the rust toolchain using `rustup`.
This will install the rust compiler `rustc` and the build system `cargo`.

This projects transitively depends on the crate  `gmp-mpfr-sys` for arbitrary length
integer support.
It is therefore required to have the `diffutils`, `gcc`, `m4`, `make` installed on Linux.
More information, as well as other operating systems can be found [here](https://docs.rs/gmp-mpfr-sys/1.4.7/gmp_mpfr_sys/index.html).

Run the following command to build the project.

```
cargo build --release
```

The resulting binary can be found in `./target/release/`.

You can start the solver with the following command using `cargo`
```
cargo run --release -- solve path/to/input-file
```
or using the binary directly:
```
./target/release/lia_equation_solver solve path/to/input-file
```

# Run

For maximal performance (during benchmarking), please add the flags `--silent` and `--no-solution`.

```
lia_equation_solver 0.1.0

USAGE:
    lia_equation_solver <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    help         Prints this message or the help of the given subcommand(s)
    solve        Solve system of equations
    to-smtlib    Convert input to smtlib format
```
