# x86 BLAS Level 1 Subset

A teaching project for implementing a subset of **BLAS Level 1** routines in **x86 assembly**.

## Overview

This project is intended for **students in assembly language courses** to practice:

- Writing **float, double, and long double** versions of basic linear algebra functions.  
- Using **SSE instructions** for single and double precision.  
- Using **x87 FPU instructions** for long double precision.  
- Understanding vectorization and CPU instruction sets through practical implementation.  

The library is implemented in `src/` and can be tested with the provided **unit testbench** in `test/`.

## Features

- BLAS Level 1 subset (vector operations)
- Multiple precision support: `float`, `double`, `long double`
- SSE and x87 FPU implementations
- Unit-tested with configurable verbosity

## Building

The project uses **Autotools**. Before building, the configure script must be executed, as `./configure`


After which, the project can be build with `make`

## Running Tests

The testbench program in `test/` can be run as:

```
./test/testbench [-v|-vv|-vvv] [function_name...]
```

* `-v` : verbose output for failed tests
* `-vv` : more detailed output (print partial vector contents)
* `-vvv` : maximum verbosity (print full vector contents)

`function_name...` may be a list of one or more functions to test:
* asumf/asum/asuml
* axpyf/axpy/axpyl
* dotf/dot/dotl
* nrm2f/nrm2/nrm2l
* scalf/scal/scall

If omitted, all tests will be executed.

