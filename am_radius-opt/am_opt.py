"""
Find methods with optimal linear radius of absolute monotonicity.
This module duplicates the functionality of RK-Opt/am_radius-opt
(but has the advantage of being written in Python).

Mathematically, it solves certain sequences of linear programming
problems.
"""
from __future__ import absolute_import
import numpy as np
import cvxpy as cvx
from six.moves import range

def ssp_explicit_lmm(r, tol, steps, order, return_coeffs=False):
    """
    Determine whether there exists an SSP explicit
    linear multistep method with SSP coefficient >= r.
    """
    k = steps
    p = order
    num_coeff = 2 * k  # Total number of method coefficients

    coeffs    = cvx.Variable(num_coeff)
    objective = cvx.Minimize(sum(coeffs))

    d = np.zeros( (p+1, 1) )
    B = np.zeros( (p+1, num_coeff) )

    #while max_radius - min_radius > precision:
    for i in range(p+1):
        d[i] = 1.
        for j in range(k):
            if i==j==0:
                B[i,j] = r/k**i
            else:
                B[i,j] = (r*j + i) * (float(j)/k)**(i-1) / k
            
            B[i, k+j] = (float(j)/k)**i


    constraints = [B * coeffs == d, coeffs >= 0]

    problem = cvx.Problem(objective, constraints)
    status = problem.solve(solver=cvx.ECOS,abstol=tol,reltol=tol,feastol=tol)#,verbose=True,max_iters=1000)
    if return_coeffs:
        beta = np.array(coeffs.value[:k]).ravel()
        alpha = np.array(coeffs.value[k:]).ravel() + r * beta
        return alpha, beta
    else: # Just return true if the problem is feasible
        return status != np.inf

def find_optimal_explicit_lmm(steps, order, precision = 1.e-15, tol = 1.e-12):
    """
    Find the optimal SSP explicit LMM with specified number of steps
    and order of accuracy.
    """
    from .utils import bisect
    from nodepy import lm

    max_radius = 1.    # upper bound on radius of a.m.

    # Find largest feasible radius
    r = bisect(0,max_radius,precision,tol,ssp_explicit_lmm,steps=steps,order=order)

    # Now find coefficients of the optimal method
    alpha, beta = ssp_explicit_lmm(r, tol, steps, order, return_coeffs=True)
    # Construct LMM
    alpha = np.concatenate((-alpha,[1]))
    beta  = np.concatenate((beta,[0]))
    return r, lm.LinearMultistepMethod(alpha,beta)


