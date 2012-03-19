********
Overview
********

`RK-opt <https://github.com/ketch/RK-opt>`_ is a MATLAB package designed to find 
Runge-Kutta (RK) methods with some optimal time stepping properties for a given stability
function :math:`R(z)`, the order of accuracy :math:`p` and its 
number of stages :math:`s` or just the last two parameters (i.e. :math:`p` and 
:math:`s`). It can compute the coefficients of explicit, 
diagonally implicit and fully implicit RK schemes. This package 
uses the MATLAB's fmincon function to handle nonlinear optimization problems. 

MATLAB global optimization toolbox with *Multistart* can be used to exploit the 
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

The following packages are now deprecated:
    * **SSP**: given the order of the scheme :math:`p` and its number of stages :math:`s` it can find the optimal RK method in terms of **SSP coefficient**
    * **low-storage**: given the order of the scheme :math:`p` and its number of stages :math:`s` it can find the optimal RK method that has the **minimal leading truncation error coefficient**

Some of the functionality in these deprecated packages has not yet been incorporated
in **RK-coeff-opt**.

A Python version of the package is also planned.


*******
RK-opt
*******

.. toctree::
   :maxdepth: 2

   started
   structure_general
   tutorial
   about
   future
