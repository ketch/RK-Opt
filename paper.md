---
title: '`RK-Opt`: A package for the design of numerical ODE solvers'
tags:
  - Python
  - numerical analysis
  - differential equations
  - Runge-Kutta method
  - linear multistep method
  - strong stability preservation
  - absolute stability
authors:
  - name: David I. Ketcheson^[Corresponding author.]
    orcid: 0000-0002-1212-126X
    affiliation: 1
  - name: Matteo Parsani
    orcid: 0000-0001-7300-1280
    affiliation: 1
  - name: Aron Ahmadia
    orcid: 0000-0002-2573-2481
    affiliation: 2
  - name: Hendrik Ranocha
    orcid: 0000-0002-3456-2277
    affiliation: 1
affiliations:
 - name: King Abdullah University of Science and Technology
   index: 1
 - name: Capital One
   index: 2
date: 9 July 2020
bibliography: paper.bib
---

# Summary

Ordinary and partial differential equations (ODEs and PDEs) are used to model
many important phenomena.  In most cases, solutions of these models must be
approximated by numerical methods.  Most of the relevant algorithms fall within
a few classes of methods, with the properties of individual methods determined
by their coefficients.  The choice of appropriate coefficients in the design of
methods for specific applications is an important area of research.
`RK-Opt` is a software package for designing numerical ODE solvers with
coefficients optimally chosen to provide desired properties.
The original (and still primary) focus of the package is on the design of
Runge-Kutta methods, but some routines for designing other classes of methods
are also included.

# Statement of need

Over the last several decades, a great deal of work has gone into the design
of numerical ODE solvers.  Initially this work was aimed at developing general
purpose solvers, but over time the emphasis shifted increasingly toward
development of optimized methods for specific applications.  Different
accuracy, stability, performance, and other properties may be relevant or
essential depending on the nature of the equations to be solved.

`RK-Opt` provides code that can enforce desired properties and/or objective
functions.  The constraints and objective are then used within an optimization
framework, to determine coefficients of methods that best achieve the desired
goal.  Thus `RK-Opt` is a sort of meta-software, consisting of algorithms whose
purpose is to create other algorithms.

Typically, the most obvious formulation of the corresponding
optimization is intractable.  Therefore, these optimization problems
are reformulated in ways that make them amenable to available techniques.
These reformulations include, for instance, turning a nonconvex problem into
a sequence of convex problems or even linear programs.  The resulting algorithms
can often guarantee optimality of their output.  However, for the
general problem of determining Runge-Kutta coefficients, the nonconvex problem
must be attacked directly and optimality cannot be guaranteed.

`RK-Opt` is written entirely in MATLAB, and leverages the MATLAB Optimization
Toolbox as well as the Global Optimization Toolbox.
Its development has been motivated largely by research needs and
it has been used in a number of papers (see below).

# Features

`RK-Opt` includes the following subpackages.

## `polyopt`

This package computes optimal stability functions for Runge-Kutta methods.
Here optimal means that the stable step size is maximized for a given ODE
spectrum.  The corresponding optimization problem is intractable under a
direct implementation.  The package uses the algorithm developed in
[@2012_optimal_stability_polynomials], which transforms the problem into a
sequence of convex problems and typically yields a solution in a few seconds or
less.  This package is usually used as the first step in designing a
Runge-Kutta method.

## `RK-Coeff-Opt`

This package computes optimal Runge-Kutta coefficients based on a desired
set of constraints and an objective.  Available constraints include:

 - The number of stages and order of accuracy
 - The class of method (explicit, implicit, diagonally implicit, low-storage)
 - The coefficients of the stability polynomial (usually determined using `polyopt`)

Two objective functions are provided; methods can be optimized for the
strong stability preserving (SSP) coefficient or the principal error norm
(a measure of the leading-order truncation error coefficients).
In addition to standard Runge-Kutta methods, various classes of multistep
Runge-Kutta methods can also be optimized.

The optimization problem in question is highly nonconvex and the available
solvers may fail to find a solution or converge to a non-optimal solution.
For this reason, the implementation is based on doing many local optimizations
in parallel from different random initial points, using MATLAB's Global
Optimization Toolbox.

The packages `dwrk-opt` and `low-storage` are specialized but less full-featured
versions of `RK-Coeff-Opt` that were developed for specific research projects
involving downwind Runge-Kutta methods and low-storage Runge-Kutta methods.

## `am_radius-opt`

Whereas the previous two subpackages are fairly general-purpose tools,
this package solves a very specific set of problems described in
[@2009_monotonicity].  Specifically, the provided routines determine the coefficients of
multistep methods (including classes of general linear methods) with the
largest possible SSP coefficient (also known
as radius of absolute monotonicity).  The corresponding optimization problem
had previously been attacked using brute force search, but this limited
its solvability to methods with very few steps.  In this package the
problem is reformulated as a sequence of linear programming problems,
enabling its efficient solution for methods with many steps.


# Related research and software

`RK-Opt` development has proceeded in close connection to the `NodePy` package (https://github.com/ketch/NodePy).
Whereas `RK-Opt` is focused on the design of numerical methods, `NodePy` is focused
more on their analysis.  A common workflow involves generating new methods with
`RK-Opt` and then studying their properties in more detail using `NodePy`.

Some of the research projects that have made use of `RK-Opt` include development of:

 - SSP Runge-Kutta methods
   [@2008_explicit_ssp;@2009_implicit_ssp;@gottlieb2015optimal]
 - SSP linear multistep methods [@2009_monotonicity]
 - SSP general linear methods [@2011_tsrk;@2017_msrk]
 - SSP IMEX Runge-Kutta methods [@conde2017implicit]
 - Low-storage Runge-Kutta methods [@2010_LSRK]
 - Optimal Runge-Kutta stability polynomials [@2012_optimal_stability_polynomials]
 - Additive and downwind SSP Runge-Kutta methods [@2011_dwssp;@2018_perturbations]
 - Optimal Runge-Kutta methods for specific PDE semi-discretizations [@parsani-eccomas;@Parsani_finnish;@2013_sd_erk;@2014_ssp_rkdg]
 - Optimal Runge-Kutta methods for pseudo-time stepping [@vermeire2019optimal;@vermeire2020optimal]
 - Embedded pairs for Runge-Kutta methods [@conde2018embedded]
 - Runge-Kutta methods with high weak stage order [@2018_wso]
 - SSP multistage, multiderivative methods [@christlieb2016explicit;@grant2019strong;@reynoso2017strong]

As can be seen from this list, applications have mostly stemmed from the
work of the main developer's research group, but have since expanded
beyond that.

# Acknowledgements

Much of the initial `RK-Opt` development was performed by D. Ketcheson while
he was supported by a DOE Computational Science Graduate Fellowship.  Development
has also been supported by funding from King Abdullah University of Science and Technology.

# References
