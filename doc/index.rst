********
Overview
********

`RK-opt <https://github.com/ketch/RK-opt>`_ is a MATLAB package for designing
Runge-Kutta (RK) methods and stability polynomials.
Supported objective functions include the principal
error norm and the SSP coefficient.  Supported constraints include stability
polynomial coefficients, low-storage formulations, and structural constraints
(explicit, diagonally implicit, etc.)
RK-opt uses MATLAB's optimization toolbox, in particular *fmincon* and *linprog*.

MATLAB's global optimization toolbox function *Multistart* can be used to exploit the
benefits of parallel search on multicore machines.


The RK-Opt package consists of the following packages:
    + :ref:`RK-coeff-opt`:
                        Find optimal Runge-Kutta method coefficients for a prescribed order of accuracy and number of stages.
                        The objective function
                        can be chosen as either the **SSP coefficient** or the
                        **leading truncation error coefficient**.
                        The method may be constrained to have a **low-storage implementation**
                        and/or a prescribed **stability polynomial**.
                        Implicit and diagonally implicit methods can also be optimized.
    + :ref:`am_radius-opt`:
                        Find stability functions with optimal radius of absolute monotonicity.
                        This includes codes for optimizing stability functions of
                        multistep, multistage methods and even methods with downwinding.
                        The optimization of rational functions is experimental.
    + :ref:`polyopt`:
                        Given a spectrum (typically corresponding to a spatial
                        semi-discretization of a PDE), find an optimal stability
                        polynomial in terms of its coefficients.  These polynomial
                        coefficients can then be used
                        as input to **RK-coeff-opt** to find a corresponding Runge-Kutta
                        method.
    + :ref:`RKtools`:
                        Some general utilities for analyzing Runge-Kutta methods.

RK-opt has been developed by David Ketcheson (primary developer and maintainer),
Matteo Parsani, and Aron Ahmadia.  It is released under a modified BSD License.
If you use RK-Opt in published work, please see :ref:`citing`. 

*******
RK-opt
*******

.. toctree::
   :maxdepth: 2

   rk-methods
   started
   citing
   applications

*********
Reference
*********
This section contains a compilation of the documentation of each function,
organized by subpackage.

.. toctree::
   :maxdepth: 2

   RK-coeff-opt
   am_radius-opt
   polyopt
   RKtools


************
Contributing
************

If you wish to contribute, we recommend that you
fork the `RK-Opt GitHub repository <https://github.com/ketch/RK-opt>`_,
implement your additions, and `issue a pull request
<https://help.github.com/articles/using-pull-requests>`_.  You may also
simply e-mail a patch to us.
