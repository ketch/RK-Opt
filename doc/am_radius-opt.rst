.. _am_radius-opt:

=============
am_radius-opt
=============
Find stability functions with optimal radius of absolute monotonicity.
This includes codes for optimizing stability functions of 
multistep, multistage methods and even methods with downwinding.

Generally, the optimization problem is phrased as a sequence of linear 
programming feasibility problems.  For details, see :cite:`2009_monotonicity`.

The optimization of rational functions is experimental.



.. contents::

multi_R_opt
===================================================
::

    function multi_R = multi_R_opt(k,p,class,varargin)


This function is a script to run the routines Rskp, Rkp_dw, Rkp_imp, or
Rkp_imp_dw several times with different inputs, in order to construct tables
of optimal values like those that appear in :cite:`2009_monotonicity`.

The inputs k, p, and (optionally) s should be vectors containing
the numbers of steps, orders of accuracy, and numbers of stages
to be considered, respectively.  The output includes results for
all combinations of values from the input vectors.

The family of methods to be considered is specified in the string 'class'.
Valid values are:

 * 'skp': find optimal general linear methods (multistep, multistage).
          In this case the vector s must be included in the inputs.
 * 'kp_imp': find optimal implicit linear multistep methods.
 * 'kp_dw':  find optimal explicit downwind linear multistep methods.
 * 'kp_imp_dw':  find optimal implicit downwind linear multistep methods.



Rkp
=================================
::

    function [R,alpha,beta]=Rkp(k,p)


Find the optimal SSP k-step explicit LMM with order of accuracy p.

Inputs: 
      * k = # of steps
      * p = order of accuracy

Outputs: 
      * R = the SSP coefficient
      * alpha, beta = the coefficients of the method

Requires MATLAB's optimization toolbox for the LP solver.



radimpfast
=============================
::

    function rad=radimpfast(p,q)


Compute the radius of absolute monotonicity of the rational function
whose numerator has coefficients p and denominator has coefficients q.
The coefficients are ordered in ascending powers.

This function is outdated and needs to be fixed.

Uses van de Griend's algorithm :cite:`vandegriend1986`, assuming multiplicity
one for all roots.  Uses high precision arithmetic.



Rskp
===============================
::

    function [R,gamma]=Rskp(s,k,p)


Finds the optimal contractive k-step, s-stage GLM with order of accuracy p
for linear problems

Inputs: s = # of stages
        k = # of steps
        p = order of accuracy

Outputs: 
       R = threshold factor
       gamma = coefficients of the polynomials
        
       for k=1, the resulting polynomial is
       `\sum_{j=0}^m (1+z/R)^j`

       For details on the general case, see :cite:`2009_monotonicity`.

This routine requires MATLAB's optimization toolbox for the LP solver.



Rkp_dw
==========================================
::

    function [R,alpha,beta,tbeta]=Rkp_dw(k,p)


Finds the optimal SSP k-step explicit LMM with order of accuracy p
allowing downwind operators

Inputs: 
      * k = # of steps
      * p = order of accuracy

Outputs: 
      * R = the SSP coefficient
      * alpha, beta, tbeta = the coefficients of the method

The method is given by
`u_n = \sum_{j=0}^{k-1} (\alpha[j] + \beta[j] F(u_{n-k+j}) + tbeta[j] tF(u_{n-k+j}))`
where tF(u) is the negated downwind operator.

Depends on MATLAB's optimization toolbox for the LP solver



Rkp_imp
=====================================
::

    function [R,alpha,beta]=Rkp_imp(k,p)


Find the optimal SSP k-step implicit LMM with order of accuracy p

Inputs: 
      * k = # of steps
      * p = order of accuracy

Outputs:
      * R = the SSP coefficient
      * alpha, beta = the coefficients of the method

Depends on MATLAB's optimization toolbox for the LP solver



Rkp_imp_dw
========================================
::

    function [R,alpha,beta]=Rkp_imp_dw(k,p)


Finds the optimal k-step implicit LMM with order of accuracy p
allowing downwinding

Inputs: 
      * k = # of steps
      * p = order of accuracy

Outputs: 
      * R = the SSP coefficient
      * alpha, beta, tbeta = the coefficients of the method
   
Depends on MATLAB's optimization toolbox for the LP solver






