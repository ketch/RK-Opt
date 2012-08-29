============
RK-coeff-opt
============

.. contents::

check_RK_order
===================================
::

    function p = check_RK_order(A,b,c)

Determines order of a RK method, up to sixth order.

For an s-stage method, input `A` should be a `s \times s` matrix;
`b` and `c` should be column vectors of length `s`.



errcoeff
===============================
::

    function D = errcoeff(A,b,c,p)


**Inputs**:
   - `A`, `b`, `c` -- Butcher tableau
   - `p`         -- order of accuracy of the method

Computes the norm of the vector of truncation error coefficients
for the terms of order `p+1`: 
(elementary weight - 1/(density of the tree)/(symmetry of the tree)


For now we just use Butcher's approach.  We could alternatively use Albrecht's.






linear_constraints
===================================================================
::

    function [Aeq,beq,lb,ub] = linear_constraints(s,class,objective,k)


This sets up:

      * The linear constraints:                     Aeq*x = beq
      * The upper and lower bounds on the unknowns: ub, lb



nonlinear_constraints
==================================================================================================
::

    function [con,coneq]=nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k)

Impose nonlinear constraints:
  - if objective = 'ssp' : both order conditions and absolute monotonicity conditions
  - if objective = 'acc' : order conditions
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
explored carefully. The Albrecht order conditions are implemented up to order 9, assuming
a certain stage order, while the Butcher order conditions are implemented up to order 9 but
do not assume anything about the stage order. Albrecht's approach is used
by default.



oc_albrecht
====================================
::

    function coneq=oc_albrecht(A,b,c,p)


Order conditions for SSP RK Methods.

This version is based on Albrecht's approach



oc_butcher
===================================
::

    function coneq=oc_butcher(A,b,c,p)


Order conditions for RKMs.
This version is based on Butcher's approach.

Assumes `p>1`.



oc_ksrk
=======================================
::

    function coneq= oc_ksrk(A,b,D,theta,p)

Order conditions for multistep-RK methods.



order_conditions
=====================================================
::

    function tau = order_conditions(x,class,s,p,Aeq,beq)


This is just a small wrapper, used when solveorderconditions=1.



rk_obj
=============================================
::

    function [r,g]=rk_obj(x,class,s,p,objective)

Objective function for RK optimization.

The meaning of the input arguments is as follow:
    * :math:`x`: vector of the unknowns.
    * class: class of method to search ('erk' = explicit RK; 'irk' = implicit RK; 'dirk' = diagonally implicit RK; 'sdirk' = singly diagonally implicit RK; '2S', '3S', '2S*', '3S*' = low-storage formulations).
    * :math:`s`:number of stages.
    * :math:`p`: order of the RK scheme.
    * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient).

The meaning of the output arguments is as follow:
    * r: it is a scalar containing the radius of absolute monotonicity if objective = 'ssp' or the value of the leading truncation error coefficient if objective = 'acc'.
    * g: it is a vector and contains the gradient of the objective function respect to the unknowns.  It is an array with all zero elements except for the last component which is equal to one if objective = 'ssp' or it is an empty array if objective = 'acc'. 



rk_opt
===================================================
::

    function rk = rk_opt(s,p,class,objective,varargin)


Find optimal RK and multistep RK methods.
The meaning of the arguments is as follows:

    * :math:`s`: number of stages.
    * :math:`k`: number of steps (1 for RK methods)
    * :math:`p`: order of the Runge-Kutta (RK) scheme.
    * class: class of method to search %       Available classes:

      * 'erk'      : Explicit Runge-Kutta methods
      * 'irk'      : Implicit Runge-Kutta methods
      * 'dirk'     : Diagonally implicit Runge-Kutta methods
      * 'sdirk'    : Singly diagonally implicit Runge-Kutta methods
      * '2S', etc. : Low-storage explicit methods; see *Ketcheson, "Runge-Kutta methods with minimum storage implementations". J. Comput. Phys. 229(5):1763 - 1773, 2010*)
      * 'emsrk1/2'    : Explicit multistep-Runge-Kutta methods
      * 'imsrk1/2'    : Implicit multistep-Runge-Kutta methods
      * 'dimsrk1/2'   : Diagonally implicit multistep-Runge-Kutta methods

    * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient)
      Accuracy optimization is not currently supported for multistep RK methods
    * poly_coeff_ind: index of the polynomial coefficients to constrain (:math:`\beta_j`) for :math:`j > p`  (j denotes the index of the stage). The default value is an empty array.  Note that one should not include any indices :math:`i \le p`, since those are determined by the order conditions.
    * poly_coeff_val: constrained values of the polynomial coefficients (:math:`\beta_j`) for :math:`j > p` (tall-tree elementary weights). The default value is an empty array.
    * startvec: vector of the initial guess ('random' = random approach; 'smart' = smart approach; alternatively, the user can provide the startvec array. By default startvec is initialize with random numbers.
    * solveorderconditions: if set to 1, solve the order conditions first before trying to optimize. The default value is 0.
    * np: number of processor to use. If np :math:`> 1` the MATLAB global optimization toolbox *Multistart* is used. The default value is 1 (just one core).
    * max_tries: maximum number of fmincon function calls. The default value is 10.
    * writeToFile: whether to write to a file. If set to 1 write the RK coefficients to a file called "ERK-p-s.txt". The default value is 1.
    * algorithm: which algorithm to use in fmincon: 'sqp','interior-point', or 'active-set'. By default sqp is used.

    .. note::
       **numerical experiments have shown that when the objective function is the minimization of the leading truncation error coefficient, the interior-point algorithm performs much better than the sqp one.**
    
    * display: level of display of fmincon solver ('off', 'iter', 'notify' or 'final'). The default value is 'notify'.
    * problem_class: class of problems for which the RK is designed ('linear' or 'nonlinear' problems). This option changes the type of order conditions check, i.e. linear or nonlinear order conditions controll. The default value is 'nonlinear'.


.. note::

   Only :math:`s`: , :math:`p`: , class and objective are required inputs.
   All the other arguments are **parameter name - value arguments to the input 
   parser scheme**. Therefore they can be specified in any order.

   **Example**::

    >>> rk=rk_opt(4,3,'erk','acc','max_tries',2,'np',1,'solveorderconditions',1)

The fmincon options are set through the **optimset** that creates/alters optimization options structure. By default the following additional options are used:
    * MaxFunEvals = 1000000
    * TolCon = 1.e-13
    * TolFun = 1.e-13
    * TolX = 1.e-13
    * MaxIter = 10000
    * Diagnostics = off
    * DerivativeCheck = off
    * GradObj = on, if the objective is set equal to 'ssp'



set_n
==========================
::

    function n=set_n(s,class)

Set total number of decision variables



shuosher2butcher
===============================================
::

    function [A,b,c]=shuosher2butcher(alpha,beta);


Generate Butcher form of a Runge-Kutta method,
given its Shu-Osher or modified Shu-Osher form.

For an m-stage method, `\alpha` and `\beta` should be 
matrices of dimension `(m+1) \times m`.



test_SSP
===============================
::

    function test_suite = test_SSP

A set of verification tests for the RK-opt package.
Currently this tests SSP coefficient optimization and
accuracy optimization, but not under constraints on the
stability polynomial.



unpack_lsrk
==============================================================================
::

    function [A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,s,class)


Extracts the coefficient arrays from the optimization vector.

This function also returns the low-storage coefficients.



unpack_msrk
=============================================================
::

    function [A,Ahat,b,bhat,D,theta] =  unpack_msrk(X,s,k,class)


Extract the coefficient arrays from the optimization vector



unpack_rk
======================================
::

    function [A,b,c]=unpack_rk(X,s,class)


Extracts the coefficient arrays from the optimization vector.

The coefficients are tored in a single vector x as::

      x=[A b' c']

A is stored row-by-row.



writeField
============================================
::

    function wf=writeField(writeFid,name,value)




