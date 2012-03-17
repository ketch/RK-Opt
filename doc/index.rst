********
Overview
********

`RK-Opt <https://github.com/ketch/RK-opt>`_ is a MATLAB package designed to find 
Runge-Kutta methods with optimal time stepping properties for a given stability
function :math:`R(z)`, the order of the scheme :math:`p` and its number of 
stages :math:`s`. This package uses the MATLAB's fmincon function to handle 
nonlinear optimization problems. 

MATLAB global optimization toolbox with *Multistart* can be used to exploit the 
benefits of parallel search on multicore machines.

RK-Opt features:
    
    * Find the optimal Runge-Kutta method with the **minimal leading truncation error coefficient**
    * Find Runge-Kutta scheme in terms of **optimal SSP coefficient**

For both the aforementioned objective functions, RK-Opt allows to impose 
low-storage properties. 

A Python version of the package is planned.


******
RK-Opt
******

.. toctree::
   :maxdepth: 2

   started
   tutorial
   about
   future
