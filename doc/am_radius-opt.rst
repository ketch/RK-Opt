=============
am_radius-opt
=============
Find stability functions with optimal radius of absolute monotonicity.
This includes codes for optimizing stability functions of 
multistep, multistage methods and even methods with downwinding.

Generally, the optimization problem is phrased as a sequence of linear 
programming feasibility problems.  For details, see [ketcheson2009]_.

The optimization of rational functions is experimental.



.. contents::




multi_R_opt
===================================================
::

    function multi_R = multi_R_opt(k,p,class,varargin)


This function is a script to run the routines Rskp, Rkp_dw, Rkp_imp, or
Rkp_imp_dw several times with different inputs, in order to construct tables
of optimal values like those that appear in [ketcheson2009]_.
different values of the input parameters, i.e.: 

k = [k1, k2, ..., kK]^T, K = length(k),  ith-element = # of steps
p = [p1, p2, ..., pP]^T, P = length(p),  ith-element = order of accuracy

and 

s = [s1, s2, ..., sS]^T, S = length(s),  ith-element = # of stages

when optimal contractive k-step, s-stage GLM are investigated.

The family of method to be considered is specified in the string 'class'.

Note that in general `S\ne K\ne P`. Fixed the order of accuracy of the time 
   integration scheme, one is usually interested in understanding the
   behavior of the threshold factor R as a function of the number of
   stages. Therefore, for a fixed element of the array "p", this function
   loops over the elements of the array "s". Thus, min(s) => max(p). The
   equality holds for any order of accuracy because the number of 
   linear order conditions that will be imposed to construct the 
   GLM coefficients is p. 



radimpfast
=============================
::

    function rad=radimpfast(p,q)


Compute the radius of absolute monotonicity of a rational function.

This function is outdated and needs to be fixed.

Uses van de Griend's algorithm, assuming multiplicity one for all roots.
Uses high precision arithmetic.



Rkp
=================================
::

    function [R,alpha,beta]=Rkp(k,p)


Find the optimal SSP k-step explicit LMM with order of accuracy p.

Inputs: 
      * `k` = # of steps
      * `p` = order of accuracy

Outputs: 
      * `\alpha, \beta` = the coefficients of the method

Requires MATLAB's optimization toolbox for the LP solver.



Rkp_dw
==========================================
::

    function [R,alpha,beta,tbeta]=Rkp_dw(k,p)


Finds the optimal SSP k-step explicit LMM with order of accuracy p
allowing downwind operators

Inputs: k = # of steps
        p = order of accuracy

Outputs: alpha, beta, tbeta = the coefficients of the method

The method is given by
`u_n = \sum_{j=0}^{k-1} (\alpha[j] + \beta[j] F(u_{n-k+j} + tbeta[j] tF(u_{n-k+j}))`
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

Outputs: alpha, beta = the coefficients of the method

Depends on MATLAB's optimization toolbox for the LP solver



Rkp_imp_dw
========================================
::

    function [R,alpha,beta]=Rkp_imp_dw(k,p)


Finds the optimal k-step implicit LMM with order of accuracy p
allowing downwinding

Inputs: k = # of steps
       p = order of accuracy

Outputs: alpha, beta, tbeta = the coefficients of the method
   
Depends on MATLAB's optimization toolbox for the LP solver



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

       in general, the resulting stability function is
       (Fill in)

epends on MATLAB's optimization toolbox for the LP solver



