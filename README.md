# Linear Integer Arithmetic Equation solver

## Build
To build this solver it is recommended to use rust toolchain using `rustup`.
This will install the rust compiler `rustc` and the build system `cargo`.

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