.. _structure_lowstorage:


===================================
RK-opt/low-storage structure solver
===================================
The main script of the **low-storage solver** is **gopt_lsrk.m** which is stored in 
**RK-opt/low-storage/**. 
Given the order of accuracy :math:`p` of the scheme and its number of stages 
:math:`s`, ssprk_opt setups and calls the MATLAB's 
`fmincon <http://www.mathworks.com/help/toolbox/optim/ug/fmincon.html>`_ 
routine to find a Runge-Kutta method that satisfies the order conditions 
(nonlinear equality constraints) and the low-storage properties 
(objective function).

gopt_lsrk MATLAB script
-----------------------
The gopt_lsrk script can be execute in MATLAB by just providing the order of the
scheme :math:`p` and its number of stages :math:`s`. 

The imposition of the objective function and nonlinear 
inequality constraints is done by calling the functions **lowstorage_obj** and 
**oc_lowstorage**.


