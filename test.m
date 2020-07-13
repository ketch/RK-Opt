% Test.m Test all sub-packages of RK-Opt

result_polyopt = runtests('polyopt/test_polyopt.m');

results_rkopt = runtests('RK-coeff-opt/test_rkopt.m');


table(result_polyopt)
table(results_rkopt)
