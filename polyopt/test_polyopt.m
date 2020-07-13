% test_polyopt.m A set of verification tests for the polyopt suite.

% setup
cvx_clear;
tol = 1.e-2;

%% Test real axis monomial
lambda = spectrum('realaxis', 100);
s = 5;
p = 1;
[h, poly_coeff] = opt_poly_bisect(lambda, s, p, 'monomial');
assert(all(ismembertol(h, 2*s^2, tol, 'ByRows', true)))

%% Test real axis Chebyshev
lambda = spectrum('realaxis', 100);
s = 8;
p = 1;
[h, poly_coeff] = opt_poly_bisect(lambda, s, p, 'chebyshev');
assert(all(ismembertol(h, 2*s^2, tol, 'ByRows', true)))

%% Test imaginary axis monomial
lambda = spectrum('imagaxis', 100);
s = 4;
p = 1;
[h, poly_coeff] = opt_poly_bisect(lambda, s, p, 'monomial');
assert(all(ismembertol(h, s-1, tol, 'ByRows', true)))

%% Test imaginary axis Chebyshev
% For some reason, SDPT3 and Sedumi sometimes fail for small s
lambda = spectrum('imagaxis', 100);
s = 8;
p = 1;
[h, poly_coeff] = opt_poly_bisect(lambda, s, p, 'rotated chebyshev')
assert(all(ismembertol(h, s-1, tol, 'ByRows', true)))

%% Test disk monomial
lambda = spectrum('disk', 100);
s = 5;
p = 1;
[h, poly_coeff] = opt_poly_bisect(lambda, s, p, 'monomial')
assert(all(ismembertol(h, s, tol, 'ByRows', true)))

%% Test disk binomial
% For some reason, SDPT3 and Sedumi sometimes fail for small s
lambda = spectrum('disk', 100);
s = 10;
p = 2;
[h, poly_coeff] = opt_poly_bisect(lambda, s, p, 'binomial')
assert(all(ismembertol(h, s-1, tol, 'ByRows', true)))
