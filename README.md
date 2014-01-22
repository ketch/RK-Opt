
# RK-opt

See the full documentation [here](http://numerics.kaust.edu.sa/RK-opt).

RK-opt is a collection of MATLAB code for designing optimized Runge-Kutta methods.
It is primarily developed and used by the 
[KAUST Numerical Mathematics Group](http://numerics.kaust.edu.sa).
It includes the following sub-packages:

 - **RK-coeff-opt**: Find optimal Runge-Kutta method coefficients, for a prescribed
   order of accuracy and number of stages.
 - **am_rad-opt**: Find stability functions with optimal radius of absolute monotonicity.
   Includes capabilities for both multistep and multistage methods.
 - **polyopt**: Find optimal stability polynomials of a given degree and order of
   accuracy for a specified spectrum.
 - **RKtools**: A collection of routines for analyzing or computing various 
   properties of Runge-Kutta methods.  For a much more extensive package along these
   lines, see [NodePy](http://numerics.kaust.edu.sa/nodepy).


# Authors
The following people have contributed code to RK-opt (listed alphabetically):

 - Aron Ahmadia: **polyopt** routines
 - David Ketcheson: Principal author and maintainer
 - Matteo Parsani: **RK-coeff-opt** routines and organization
 - Zachary Grant: Extension of order conditions to multistep RK with more than two stages and 
    addition of order conditions for orders 9-11.
