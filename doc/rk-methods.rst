.. _rk:

========================================
Automated design of Runge-Kutta methods
========================================
An `s`-stage Runge-Kutta method has roughly `s^2` coefficients (roughly `s^2/2` for explicit methods), 
which can be chosen
so as to provide high accuracy, stability, or other properties.  Historically, most
interest in Runge-Kutta methods has focused on methods using the minimum number of stages
for a given order of accuracy.  However, in the past few decades there has been increasing
recognition that using *extra* stages can be worthwhile in order to improve other method
properties.  Some areas where this is particularly useful are in the enhancement of linear
and nonlinear stability properties, the reduction of storage requirements, and the design of embedded pairs.
Methods with dozens or even hundreds of stages are not unheard of.

At the same time, most existing Runge-Kutta methods have been designed by hand, by researchers laboriously
solving the order conditions.  When using extra stages, the number of available parameters makes the
selection of a near-optimal choice by hand impossible, and one resorts to computational optimization.
This leads to a different paradigm of numerical method design, in which we use sophisticated numerical
(optimization) algorithms to design sophisticated numerical (integration) algorithms.  It can be expected
that this trend will accelerate in the future, and perhaps one day simple manually-constructed
algorithms will be the exception.

RK-Opt contains a set of tools for designing Runge-Kutta methods in this paradigm.  It has been
constructed mostly in the direct line of our research, but we have made some effort to help others
easily understand and use it.  We hope that you find it useful, and that you will contribute any
enhancements you may develop back to the project by sending us a pull request on Github.
