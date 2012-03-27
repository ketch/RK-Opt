.. _structure_general:


================================
RK-opt/general structure solver
================================
The main function of the **general solver** is **rk_opt** which is stored in 
**RK-opt/general/rk_opt.m**. 
Given the order of 
accuracy :math:`p` of the scheme, its number of stages :math:`s`, 
rk_opt setups and calls the MATLAB's 
`fmincon <http://www.mathworks.com/help/toolbox/optim/ug/fmincon.html>`_ 
routine to find a Runge-Kutta method that satisfies the order conditions 
(nonlinear equality constraints) and 
a combination of inequality constraints and objective functions. The method may 
be also constrained to have a **low-storage implementation** and/or a prescribed 
**stability polynomial**.

In this page, the rk_opt function and few other key routines are briefly 
described.


rk_opt function
---------------
The rk_opt function signature is ::

     rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,np,max_tries,writeToFile,algorithm,display).

The meaning of the arguments is as follow:
    * :math:`s`: number of stages.
    * :math:`p`: order of the Runge-Kutta (RK) scheme.
    * class: class of method to search ('erk' = explicit RK; 'irk' = implicit RK; 'dirk' = diagonally implicit RK; 'sdirk' = singly diagonally implicit RK; '2S', '3S', '2S*', '3S*' = low-storage formulations, see *Ketcheson, "Runge-Kutta methods with minimum storage implementations". J. Comput. Phys. 229(5):1763 - 1773, 2010*)
    * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient)
    * poly_coeff_ind: index of the polynomial coefficients (:math:`\beta_j`) for :math:`j > p`  (j denotes the index of the ). The default value is an empty array.
    * poly_coeff_val: values of the polynomial coefficients (:math:`\beta_j`) for :math:`j > p` (tall-tree elementary weights). The default value is an empty array.
    * startvec: vector of the initial guess ('random' = random approach; 'smart' = smart approach; alternatively, the user can provide the startvec array. By default startvec is initialize with random numbers.
    * solveorderconditions: if set to 1, solve the order conditions first before trying to optimize. The default value is 0.
    * np: number of processor to use. If np :math:`> 1` the MATLAB global optimization toolbox *Multistart* is used. The default value is 1 (just one core).
    * max_tries: maximum number of fmincon function calls. The default value is 10.
    * writeToFile: whether to write to a file. If set to 1 write the RK coefficients to a file called "ERK-p-s.txt". The default value is 1.
    * algorithm: which algorithm to use in fmincon (sqp or interior-point). By default sqp is used.
    * display: level of display of fmincon solver ('off', 'iter', 'notify' or 'final'). The default value is 'notify'.
    * problem_class: class of problems for which the RK is designed ('linear' or 'nonlinear' problems). This option changes the type of order conditions check, i.e. linear or nonlinear order conditions controll. The default value is 'nonlinear'.


.. note::

   Only :math:`s`: , :math:`p`: , class and objective are required inputs.
   All the other arguments are **parameter name - value arguments to the input 
   parser scheme**. Therefore they can be specified in a free order.
   Example::

    >>>> rk=rk_opt(4,3,'erk','acc','max_tries',2,'np',1,'solveorderconditions',1)


The fmincon options are set through the **optimset** that creates/alters optimization options structure. By default the following additional options are used:
    * MaxFunEvals = 1000000
    * TolCon = 1.e-13
    * TolFun = 1.e-13
    * TolX = 1.e-13
    * MaxIter = 10000
    * Diagnostics = off
    * DerivativeCheck = off
    * GradObj = on, if the objective is set equal to 'ssp'


The rk_opt function is also in charge to setup the objective function and both 
linear equality and nonlinear inequality constraints. This is achieved by 
calling the functions **rk_obj**, **linear_constraints** and **nlc** which are 
described below.


rk_obj
------
The rk_obj function signature is ::
    
    [r,g]=rk_obj(x,class,s,p,objective)

The meaning of the input arguments is as follow:
    * :math:`x`: vector of the unknowns.
    * class: class of method to search ('erk' = explicit RK; 'irk' = implicit RK; 'dirk' = diagonally implicit RK; 'sdirk' = singly diagonally implicit RK; '2S', '3S', '2S*', '3S*' = low-storage formulations).
    * :math:`s`:number of stages.
    * :math:`p`: order of the RK scheme.
    * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient).

The meaning of the output arguments is as follow:
    * r: it is a scalar containing the radius of absolute monotonicity if objective = 'ssp' or the value of the leading truncation error coefficient if objective = 'acc'.
    * g: it is a vector and contains the gradient of the objective function respect to the unknowns.  It is an array with all zero elements except for the last component which is equal to one if objective = 'ssp' or it is an empty array if objective = 'acc'. 


linear_constraints function
---------------------------
The linear_constraints function signature is ::
    
    [Aeq,beq,lb,ub] = linear_constraints(s,class,objective)

The meaning of the input arguments is as follow:
    * :math:`s`: number of stages.
    * class: class of method to search ('erk' = explicit RK; 'irk' = implicit RK; 'dirk' = diagonally implicit RK; 'sdirk' = singly diagonally implicit RK; '2S', '3S', '2S*', '3S*' = low-storage formulations).
    * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient).

The meaning of the output arguments is as follow:
    * Aeq, beq: linear constraints Aeq*x = beq. These constraints depends on the class of RK scheme.
    * lb, ub: lower and upper bounds of the unknowns (i.e. the Butcher's tableau coefficients).



nlc function
------------
The **nlc** function signature is ::

    [con,coneq]=nlc(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val)

The meaning of the input arguments is as follow:
    * :math:`x`: vector of the unknowns.
    * class: class of method to search ('erk' = explicit RK; 'irk' = implicit RK; 'dirk' = diagonally implicit RK; 'sdirk' = singly diagonally implicit RK; '2S', '3S', '2S*', '3S*' = low-storage formulations).
    * :math:`s`:number of stages.
    * :math:`p`: order of the RK scheme.
    * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient).
    * poly_coeff_ind: index of the polynomial coefficients (:math:`\beta_j`) for :math:`j > p`.
    * poly_coeff_val: values of the polynomial coefficients (:math:`\beta_j`) for :math:`j > p` (tall-tree elementary weights).

The meaning of the output arguments is as follow:
    * con: inequality constraints, i.e. absolute monotonicity conditions if objective = 'ssp' or nothing if objective = 'acc'
    * coneq: order conditions plus stability function coefficients constraints (tall-tree elementary weights)

Two forms of the order conditions are implemented: one based on **Butcher's 
approach**, and one based on **Albrecht's approach**. One or the other may lead 
to a more tractable optimization problem in some cases, but this has not been 
explored carefully. The Albrecht order conditions go up to order 9, but assume 
a certain stage order, while the Butcher order conditions go up to order 6 but
do not assume anything about the stage order. The Albrecht's approach is used
by default.


unpack_x routines
-----------------
Because fmincon requires that all decision variables be packed into a single 
vector, one of the most tedious parts of the package is keeping track of the 
correspondence between that vector and the variables of interest. This is 
implemented in the **unpack_x** routines. 

Currently two unpack_x routines are available: unpack_lsrk (lsrk = low-storage
RK) and unpack_rk. The first one computes both low-storage formulation prescribed
in class ('2S', '3S', '2S*', '3S*) and the
Butcher's tableau; the second one just calculates the Butcher's tableau.





   
