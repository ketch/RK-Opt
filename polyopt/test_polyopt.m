function test_suite = test_polyopt
%function test_suite = test_polyopt
%
% A set of verification tests for the polyopt suite.

initTestSuite;

function tol = setup
    cvx_clear;
    tol = 1.e-2;

function test_realaxis_monomial(tol)
    lamda = spectrum('realaxis',100);
    s=randi(6); p=1;
    [h,poly_coeff]=opt_poly_bisect(lamda,s,p)
    assertElementsAlmostEqual(h,2*s^2,'relative',tol);

function test_realaxis_chebyshev(tol)
    lamda = spectrum('realaxis',100);
    s=randi(10); p=1;
    [h,poly_coeff]=opt_poly_bisect(lamda,s,p,'chebyshev')
    assertElementsAlmostEqual(h,2*s^2,'relative',tol);

function test_imagaxis_monomial(tol)
    lamda = spectrum('imagaxis',100);
    s=randi(5); p=1;
    [h,poly_coeff]=opt_poly_bisect(lamda,s,p,'monomial')
    assertElementsAlmostEqual(h,s-1,'relative',tol);

function test_imagaxis_chebyshev(tol)
    % For some reason, SDPT3 and Sedumi sometimes fail for small s
    lamda = spectrum('imagaxis',100);
    s=randi(10)+10; p=1;
    [h,poly_coeff]=opt_poly_bisect(lamda,s,p,'rotated chebyshev')
    assertElementsAlmostEqual(h,s-1,'relative',tol);
    
function test_disk_monomial(tol)
    lamda = spectrum('disk',100);
    s=randi(5); p=1;
    [h,poly_coeff]=opt_poly_bisect(lamda,s,p,'monomial')
    assertElementsAlmostEqual(h,s,'relative',tol);

function test_disk_binomial(tol)
    % For some reason, SDPT3 and Sedumi sometimes fail for small s
    lamda = spectrum('disk',100);
    s=randi(10)+10; p=2;
    [h,poly_coeff]=opt_poly_bisect(lamda,s,p,'binomial')
    assertElementsAlmostEqual(h,s-1,'relative',tol);
