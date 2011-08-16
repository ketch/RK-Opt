
******
RK-opt
******

RK-opt is a package designed to search for optimized Runge-Kutta methods.
The package makes use of MATLAB's `fmincon` function, which handles nonlinear
optimization problems with equality and inequality constraints.  
Equality constraints can be imposed to enforce 

    * order conditions
    * low-storage properties
    * stability function coefficients (tall-tree elementary weights)
      
A combination of inequality constraints and certain objective functions can be
used to find methods

    * That are optimal in terms of SSP coefficient
    * That have minimal leading error coefficients

Note that these are both nonlinear properties of the RK method.
Since linear properties of the method can be optimized by dealing
directly with the stability coefficients, a good approach is to
optimize them first and then enforce optimality through equality
constraints on the tall-tree elementary weights.

Because `fmincon` requires that all decision variables be packed
into a single vector, one of the most tedious parts of the package
is keeping track of the correspondence between that vector and the
variables of interest.  This is implemented in the unpack_x routines.

Optionally, the MATLAB global optimization toolbox can be used to 
take advantage of multicore machines.

Two forms of the order conditions are implemented: one based on Butcher's
approach, and one based on Albrecht's approach.  One or the other may lead to
a more tractable optimization problem in some cases, but this has not been
explored carefully.  The Albrecht order conditions go up to order 9, but
assume a certain stage order, while the Butcher order conditions go up
to order 6 but don't assume anything about the stage order.

A Python-based rewrite of the package is planned.
