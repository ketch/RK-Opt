
# RK-Opt: A Package for the Design of Numerical ODE solvers

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rk-opt.readthedocs.io/en/latest/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-success.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4138076.svg)](https://doi.org/10.5281/zenodo.4138076)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02514/status.svg)](https://doi.org/10.21105/joss.02514)

See the full documentation [here](https://rk-opt.readthedocs.io/en/latest/).

RK-Opt is a collection of MATLAB codes for designing optimized numerical ODE solvers.
The main emphasis is on Runge-Kutta methods, but some routines deal with other classes of methods.
It is primarily developed and used by the
[KAUST Numerical Mathematics Group](http://numerics.kaust.edu.sa).
It includes the following sub-packages:

 - **polyopt**: Find optimal stability polynomials of a given degree and order of
   accuracy for a specified spectrum.
 - **RK-coeff-opt**: Find optimal Runge-Kutta method coefficients, for a prescribed
   order of accuracy and number of stages.
 - **am_rad-opt**: Find stability functions with optimal radius of absolute monotonicity.
   Includes capabilities for both multistep and multistage methods.
 - **RKtools**: A collection of routines for analyzing or computing various
   properties of Runge-Kutta methods.  For a much more extensive package along these
   lines, see [NodePy](http://nodepy.readthedocs.io/en/latest/).

A common workflow for designing Runge-Kutta methods is to use **polyopt** to find an
appropriate stability function and then **RK-coeff-opt** to determine the Runge-Kutta
method coefficients.

To run the tests, execute the MATLAB script `test.m`. This requires a relatively recent
version of MATLAB (tested with R2018a and later) with the following toolboxes.
 - MATLAB Optimization Toolbox
 - MATLAB Global Optimization Toolbox
 - CVX (http://cvxr.com/cvx/)
 - MATLAB Parallel Computing Toolbox (optional; allows faster searching for optimal methods in RK-Coeff-Opt)


# Authors
The code is primarily developed and maintained by David Ketcheson.
The following people have also made important contributions to RK-Opt (listed alphabetically):

 - Aron Ahmadia: Co-developer of **polyopt** algorithm and routines.
 - Zachary Grant: Extension of order conditions to multistep RK with more than two stages and
    addition of order conditions for orders 9-11.
 - Matteo Parsani: Many improvements to **RK-coeff-opt** routines and organization.
 - Hendrik Ranocha: General improvements and updates, including updating the test routines.
