# Binding states for deuteron
Code for exactly solving neuteron's binding states in different ways.

## Prerequisites
A fortran compiler (e.g., gfortran) is installed.

## Compile and execution
### direct/
This code directly solves the equation.
```bash
gfortran direct_solution.f90
./a.out
```

### ho_basis/
This code solves the equation by expanding the wave function with eigen functions of harmonic oscillator.
```bash
gfortran ho_basis.f90
./a.out
```
