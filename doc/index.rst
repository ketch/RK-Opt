********
Overview
********

`RK-opt <https://github.com/ketch/RK-opt>`_ is a MATLAB package for designing
Runge-Kutta (RK) methods and stability polynomials.  
Supported objective functions include the principal
error norm and the SSP coefficient.  Supported constraints include stability
polynomial coefficients, low-storage formulations, and structural constraints
(explicit, diagonally implicit, etc.)
RK-opt uses MATLAB's optimization toolbox, in particular `fmincon` and `linprog`.

MATLAB's global optimization toolbox function `Multistart` can be used to exploit the 
benefits of parallel search on multicore machines.


The RK-Opt pacakge consists of the following packages:
    * **RK-coeff-opt**: Find optimal Runge-Kutta method coefficients, for a prescribed order of accuracy and number of stages.  
                        The objective function
                        can be chosen as either the **SSP coefficient** or the
                        **leading truncation error coefficient**.
                        The method may be constrained to have a **low-storage implementation**
                        and/or a prescribed **stability polynomial**.
                        Implicit and diagonally implicit methods can also be optimized.
    * **am_rad-opt**: Find stability functions with optimal radius of absolute monotonicity.
                        This includes codes for optimizing stability functions of 
                        multistep, multistage methods and even methods with downwinding.
                        The optimization of rational functions is experimental.
    * **polyopt**:    Given a spectrum (typically corresponding to a spatial
                        semi-discretization of a PDE), finds an optimal stability
                        polynomial.  The polynomial coefficients can then be used
                        as input to **RK-coeff-opt** to find a corresponding Runge-Kutta
                        method.
    * **RKtools**:    Some general utilities for analyzing Runge-Kutta methods.
                        


*******
RK-opt
*******

.. toctree::
   :maxdepth: 2

   started
   tutorial
   about
   future

*********
Reference
*********

.. toctree::
   :maxdepth: 2

   RK-coeff-opt
   am_radius-opt
   polyopt
   RKtools

Contributing
-------------
If you wish to contribute, we recommend that you
fork the RK-Opt Github repository, implement your additions, and issue
a pull request.  You may also simply e-mail a patch to us.


