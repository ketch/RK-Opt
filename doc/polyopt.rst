=======
polyopt
=======

.. contents::













opt_poly_bisect
==================================================================
::

    function [h,poly_coeff] = opt_poly_bisect(lam,s,p,basis,varargin)


Finds an optimally stable polynomial of degree s and order p for the spectrum
lam in the interval (h_min,h_max) to precision eps.



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



