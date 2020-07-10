
# RK-opt: A Package for the Design of Numerical ODE solvers

See the full documentation [here](http://numerics.kaust.edu.sa/RK-opt).

RK-opt is a collection of MATLAB codes for designing optimized numerical ODE solvers.
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
   lines, see [NodePy](http://numerics.kaust.edu.sa/nodepy).

A common workflow for designing Runge-Kutta methods is to use **polyopt** to find an
appropriate stability function and then **RK-coeff-opt** to determine the Runge-Kutta
method coefficients.

# Authors
The code is primarily developed and maintained by David Ketcheson.
The following people have also made important constributions to RK-opt (listed alphabetically):

 - Aron Ahmadia: Co-developer of **polyopt** algorithm and routines.
 - Zachary Grant: Extension of order conditions to multistep RK with more than two stages and
    addition of order conditions for orders 9-11.
 - Matteo Parsani: Many improvements to **RK-coeff-opt** routines and organization.
 - Hendrik Ranocha: General improvements and updates, including updating the test routines.
