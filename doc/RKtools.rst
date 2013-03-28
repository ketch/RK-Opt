=======
RKtools
=======
Some general utilities for analyzing Runge-Kutta methods.



.. contents::

am_radius
=======================================
::

    function r = am_radius(A,b,c,eps,rmax)


Evaluates the Radius of absolute monotonicity
of a Runge-Kutta method, given the Butcher array.

For an `m`-stage method, `A` should be an `m \times m` matrix
and `b` should be a column vector of length m.

Accuracy can be changed by modifying the value of eps (default `10^-{10}`)
Methods with very large radii of a.m. (>50) will require
the default value of rmax to be increased.

The radius of absolute monotonicity is the largest value of `r`
such that

.. raw:: latex

   \begin{eqnarray}
   K(I+rA)^{-1} \ge & 0    \\
   rK(I+rA)^{-1}e_m \le & e_{m+1}
   \end{eqnarray}

   where $$ K = \left(\begin{array}{c} A \\ b^T \end{array}\right) $$



internal_stab_explicit_butcher
====================================================================================
::

    function [stability] = internal_stab_explicit_butcher(A,b,c,spectrum,one_step_dt,p)



This function computes and plots both intermediate and one-step internal 
stability vector of an explicit Runge-Kutta scheme given its Butcher 
tableau.

Note that for an explicit Runge-Kutta scheme the stability functions are
polynomials in the complex variable z.

Construct the intermediate stability functions \psi_j (where j is the 
index of the stage).

Note that for an explicit scheme the intermediate stability polynomial 
associated to the first stage is always 1, i.e. \psi_1 = 1.
Therefore we just compute and plot the remaining (s-1) intermediate
stability polynomials plus the one-step stability polynomial of the
Runge-Kuatta method.



L2_timestep_poly
=================================================
::

    function c = L2_timestep_poly(sdisc,p,q,eps,tol)


Find the absolutely timestep for a given combination of
linear spatial discretization and stability function.

Also (optionally) plots the region of absolute stability and the eigenvalues.

The timestep is determined to within accuracy eps (default 10^-4).

The spectral stability condition is checked to within tol (default 10^-13).



optimal_shuosher_form
=======================================================
::

    function [v,alpha,beta] = optimal_shuosher_form(A,b,c)




plotstabreg(rk,plotbounds,ls,lw)
==========================================
::

    function plotstabreg(rk,plotbounds,ls,lw)


Plots the absolute stability region
of a Runge-Kutta method, given the Butcher array



plotstabreg_func
======================================================
::

    function [dummy] = plotstabreg_func(p,q,bounds,ls,lw)


plot the absolute stability region of a one-step method,
given the stability function

Inputs:
      * p: coefficients of the numerator   of the stability function
      * q: coefficients of the denominator of the stability function 

 if q is omitted, it is assumed that the function is a polynomial
Remaining inputs are optional:
      * bounds: bounds for region to compute and plot (default [-9 1 -5 5])
      * ls:   line style (default '-r')
      * lw:   line width (default 2)



rk_stabfun
================================
::

    function [p,q] = rk_stabfun(rk)


Outputs the stability function of a RK method.
The Butcher coefficients should be stored in rk.A, rk.b.

p contains the coefficients of the numerator

q contains the coefficients of the denominator


.. raw:: latex

    $$\phi(z)=\frac{\sum_j p_j z^j}{\sum_j q_j z^j} = \frac{\det(I-z(A+eb^T))}{\det(I-zA)}.$$



semispectrum
======================================================
::

    function L = semispectrum(method,order,doplot,nx,cfl)

Plot spectra of various semi-discretizations of the advection equation

Current choices for method:
      - 'fourier':   Fourier   spectral method
      - 'chebyshev': Chebyshev spectral method
      - 'updiff':    Upwind difference operators (linearized WENO)
      - 'DG':        Discontinuous Galerkin method

The value of order matters only for the 'updiff' and 'DG' methods
and selects the order of accuracy in those cases.



