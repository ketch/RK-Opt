********
Overview
********

`RK-Opt` is a software package for designing numerical ODE solvers with
coefficients optimally chosen to provide desired properties.
It is available from https://github.com/ketch/RK-Opt, with documentation
at http://rk-opt.readthedocs.io/en/latest/.
The primary focus of the package is on the design of Runge-Kutta methods
(including both stability polynomials and full Butcher tableaus), but
some routines for designing other classes of methods such as multistep
Runge-Kutta and general linear methods are also included.
Supported objective functions include the principal
error norm and the SSP coefficient.  Supported constraints include stability
polynomial coefficients, low-storage formulations, and structural constraints
(explicit, diagonally implicit, etc.)
RK-Opt uses `CVX <http://cvxr.com/cvx/>`_ as well as MATLAB's Optimization Toolbox and Global
Optimization Toolbox.

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

RK-Opt has been developed by David Ketcheson (primary developer and maintainer),
Matteo Parsani, Aron Ahmadia, Zack Grant, and Hendrik Ranocha.

RK-Opt is released under a modified BSD License.
If you use RK-Opt in published work, please cite it; see :ref:`citing`.

*******
RK-Opt
*******

.. toctree::
   :maxdepth: 2

   rk-methods
   started
   citing
   zreferences

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
fork the `RK-Opt GitHub repository <https://github.com/ketch/RK-Opt>`_,
implement your additions, and `issue a pull request
<https://help.github.com/articles/using-pull-requests>`_.  You may also
simply e-mail a patch to us.

