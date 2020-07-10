Given a spectrum (typically corresponding to a spatial
semi-discretization of a PDE), finds an optimal stability polynomial. The
polynomial coefficients can then be used as input to `RK-coeff-opt` to find a
corresponding Runge-Kutta method.

This is the implementation of the algorithm described in [ketcheson-ahmadia]_.
The code was written by Aron Ahmadia and David Ketcheson.

To run the tests, execute the MATLAB commands
```
results_polyopt = runtests('test_polyopt.m');
table(results_polyopt)
```
