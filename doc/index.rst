********
Overview
********

`RK-Opt <https://github.com/ketch/RK-opt>`_ is a MATLAB package designed to find 
Runge-Kutta (RK) methods with some optimal time stepping properties for a given stability
function :math:`R(z)`, the order of accuracy :math:`p` and its 
number of stages :math:`s` or just the last two parameters (i.e. :math:`p` and 
:math:`s`). It can compute the coefficients of explicit, 
diagonally implicit and fully implicit RK schemes. This package 
uses the MATLAB's fmincon function to handle nonlinear optimization problems. 

MATLAB global optimization toolbox with *Multistart* can be used to exploit the 
benefits of parallel search on multicore machines.


The RK-Opt pacakge consists of three main solvers stored in three different directories:
    * **general**: given a stability function :math:`R(z)`, the order of the scheme :math:`p` and its number of stages :math:`s` it can find the optimal RK method in terms of **SSP coefficient** or that has **minimal leading truncation error coefficient**. RK-Opt allows to impose **low-storage properties** for explicit RK schemes and both objective functions.
    * **SSP**: given the order of the scheme :math:`p` and its number of stages :math:`s` it can find the optimal RK method in terms of **SSP coefficient**
    * **low-storage**: given the order of the scheme :math:`p` and its number of stages :math:`s` it can find the optimal RK method that has the **minimal leading truncation error coefficient**

Currently the solvers are completely separate. However, we will probably merge 
them in a single entity in the near future.

A Python version of the package is also planned.


*******
RK-opt
*******

.. toctree::
   :maxdepth: 2

   started
   structure_general
   structure_SSP
   structure_lowstorage
   tutorial
   about
   future
