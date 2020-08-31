Some general utilities for analyzing Runge-Kutta methods.

Some of the routines expect as input a structured array `rk`.
This structure must have the fields `A, b, c`, containing its
Butcher coefficients.  Optionally, it may represent an additive
Runge-Kutta method or an embedded pair in which case it should also have
`Ahat`, `bhat`, `chat` containing the coefficients of the secondary
method.
