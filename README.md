# Lattice Boltzmann Flow Simulation 
## Description
This was a final project for a computational physics class. It creates an animation of a simple two-dimensional 
object in a fluid flow, visualizing the speed of the fluid. This project uses fortran wrapped in python, allowing 
for performance with the ease of visualization in python. 

![sphere-flow](https://github.com/skrawczuk/lattice-boltzmann-flow-sim/blob/master/sphere.gif)

## Running
This project uses python packages ```numpy```, ```pandas```, and ```matplotlib```.
To compile the Fortran module, use:
```
! f2py -c -m lbm_step lbm_step.f90 --f90flags="-fopenmp" -lgomp
```
If you don't want to use the parallel portion, simply remove the flags: 
```
! f2py -c -m lbm_step lbm_step.f90 
```
