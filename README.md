# FORTRAN Linear System Solver

This repository contains a simple Fortran solver for linear systems of the form `Ax = b`.
It uses LU factorization with partial pivoting.

## Files

- `linear_system_toolbox.f90` - solver module
- `lse_test.f90` - small test program

## Build

With `gfortran`:

```bash
gfortran linear_system_toolbox.f90 lse_test.f90 -o lse_test
```

## Run

```bash
./lse_test
```
