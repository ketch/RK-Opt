.. _structure_SSP:


================================
RK-opt/SSP structure solver
================================
The main script of the **SSP solver** is **ssprk_opt.m** which is stored in 
**RK-opt/SSP/**. 
Given the order of accuracy :math:`p` of the scheme and its number of stages 
:math:`s`, ssprk_opt setups and calls the MATLAB's 
`fmincon <http://www.mathworks.com/help/toolbox/optim/ug/fmincon.html>`_ 
routine to find a Runge-Kutta method that satisfies the order conditions 
(nonlinear equality constraints), the absolute monotonicity conditions 
(nonlinear inequality constraints) and maximises the radius of absolute 
monotonicity :math:`r` (objective function).

For this solver, the function that uses the MATLAB global optimization toolbox 
with *Multistart* for parallel search on multicore machines is not implemented.

ssprk_opt MATLAB script
-----------------------
The ssrk_opt script can be execute in MATLAB by just providing the order of the
scheme :math:`p` and its number of stages :math:`s`. 

The imposition of the objective function and both linear equality and nonlinear 
inequality constraints is done by calling the functions **rk_am_obj**, 
**linear_constraints** and **nlc**. These functions are like those of 
the general solver (see :ref:`structure_general`) when the objective function is 
set to maximize SSP coefficient :math:`r`



