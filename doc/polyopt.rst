=======
polyopt
=======
Given a spectrum (typically corresponding to a spatial
semi-discretization of a PDE), finds an optimal stability polynomial. The
polynomial coefficients can then be used as input to RK-coeff-opt to find a
corresponding Runge-Kutta method.

This is the implementation of the algorithm described in [ketcheson-ahmadia]_.
The code was written by Aron Ahmadia and David Ketcheson.



.. contents::













opt_poly_bisect
==================================================================
::

    function [h,poly_coeff] = opt_poly_bisect(lam,s,p,basis,varargin)


Finds an optimally stable polynomial of degree s and order p for the spectrum
lam in the interval (h_min,h_max) to precision eps.

Optional arguments:

      lam_func: 
                A function used to generate the appropriate spectrum
                at each bisection step, instead of using a fixed (scaled) spectrum.
                Used for instance to find the longest rectangle of a fixed height
                (see Figure 10 of the CAMCoS paper).

Examples:

      - To find negative real axis inclusion::

              lam = spectrum('realaxis',500);       
              s = 10; p = 2;
              [h,poly_coeff] = opt_poly_bisect(lam,s,p,'chebyshev')    

      - To reproduce figure 10 of [ketcheson-ahmadia]_ ::

              lam_func = @(kappa) spectrum('rectangle',100,kappa,10)
              [h,poly_coeff] = opt_poly_bisect(lam,20,1,'chebyshev','lam_func',lam_func)
              plotstabreg_func(poly_coeff,[1])



spectrum
=============================================
::

    function lamda = spectrum(name,N,kappa,beta)


Return N discretely sampled values from certain sets in the complex plane.

Acceptable values for name:
      * 'realaxis':     `[-1,0]`
      * 'imagaxis':     `[-i,i]`
      * 'disk':         `{z : |z+1|=1}`
      * 'rectangle':    `{x+iy : -\beta \le y \le \beta, -\kappa \le x \le 0}`
      * 'Niegemann-ellipse' and 'Niegemann-circle':  See Niegemann 2011
      * 'gap':          Spectrum with a gap; see Ketcheson & Ahmadia 2012

kappa and beta are used only if name == 'rectangle'



test_polyopt
===================================
::

    function test_suite = test_polyopt


A set of verification tests for the polyopt suite.



