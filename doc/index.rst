********
Overview
********

`RK-Opt <https://github.com/ketch/RK-opt>`_ is a MATLAB package designed to find 
Runge-Kutta methods with optimal time stepping properties for a given stability
function :math:`R(z)`, the order of the scheme :math:`p` and its number of 
stages :math:`s`. It can handle explicit, diagonally implicit and fully implicit
Runge-Kutta schemes. This package uses the MATLAB's fmincon function to handle 
nonlinear optimization problems. 

MATLAB global optimization toolbox with *Multistart* can be used to exploit the 
benefits of parallel search on multicore machines.

RK-Opt. features:
    
    * find the optimal Runge-Kutta method with the **minimal leading truncation error coefficient**;
    * find Runge-Kutta scheme in terms of **optimal SSP coefficient**.

RK-Opt allows to impose **low-storage properties**, for the explicit Runge-Kutta 
class and both the aforementioned objective functions.

A Python version of the package is planned.


*******
RK-Opt.
*******

.. toctree::
   :maxdepth: 2

   started
   structure
   tutorial
   about
   future
