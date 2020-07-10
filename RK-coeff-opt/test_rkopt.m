% test_rkopt.m A set of verification tests for the RK-opt package.
%
% Currently this tests SSP coefficient optimization and
% accuracy optimization, but not under constraints on the
% stability polynomial.


% Check the solution obtained with fmincon nested in a while loop
% =============================================================

%% Test SSP32
tol = 1.e-14;
rk = rk_opt(3,2,'erk','ssp','startvec','smart','write_to_file',0);
A = [0 0 0; 0.5 0 0; 0.5 0.5 0];
r = 2.;
assert(all(ismembertol(rk.r(:)', r(:)', tol, 'ByRows', true)))
assert(all(ismembertol(rk.A(:)', A(:)', tol, 'ByRows', true)))


%% Test SSP43
rk = rk_opt(4,3,'erk','ssp','startvec','smart','write_to_file',0);
b = [1./6 1./6 1./6 1./2]';
r = 2.;
assert(all(ismembertol(rk.r(:)', r(:)', 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', 'ByRows', true)))


%% Test SDIRK SSP22
rk = rk_opt(2,2,'sdirk','ssp','startvec','random','write_to_file',0);
b = [1./2 1./2]';
r = 4.;
assert(all(ismembertol(rk.r(:)', r(:)', 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', 'ByRows', true)))


%% Test DIRK SSP22
tol = 1.e-7;
rk = rk_opt(2,2,'dirk','ssp','startvec','smart','write_to_file',0);
b = [1./2 1./2]';
r = 4.;
if abs(rk.r-r)>tol || max(max(abs(rk.b-b)))>tol
    error('Failed to find optimal DIRK SSP(2,2) method. This occasionally happens; try running the tests again.')
end
assert(all(ismembertol(rk.r(:)', r(:)', tol, 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', tol, 'ByRows', true)))


%% Test RK22 acc
tol = 1.0e-5;
x = [0.117831902493812 0.187227391881097   0.812772608118903  -0.099952256513284 -0.010000000000000];
rk = rk_opt(2,2,'erk','acc','startvec',x,'write_to_file',0);
b = [0.25 0.75]';
assert(all(ismembertol(rk.errcoeff', 1./6, 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', tol, 'ByRows', true)))


%% Test SSPTSRK2
tol = 1.e-14;
s = randi(3)+1;
rk = rk_opt(s,2,'emsrk2','ssp','k',2,'write_to_file',0,'min_amrad',s-1);
r = sqrt(s*(s-1));
assert(all(ismembertol(rk.r(:)', r(:)', tol, 'ByRows', true)))



% Check the solution obtained with fmincon called by multistart solver
% both in serial and in parallel.
% ====================================================================

%% Test SSP32 multistart
tol = 1.e-14;
rk = rk_opt(3,2,'erk','ssp','startvec','smart','write_to_file',0);
A = [0 0 0; 0.5 0 0; 0.5 0.5 0];
r = 2.;
assert(all(ismembertol(rk.r(:)', r(:)', tol, 'ByRows', true)))
assert(all(ismembertol(rk.A(:)', A(:)', tol, 'ByRows', true)))


%% Test SSP32 multistart parallel
tol = 1.e-14;
rk = rk_opt(3,2,'erk','ssp','startvec','smart','np',2,'write_to_file',0);
A = [0 0 0; 0.5 0 0; 0.5 0.5 0];
r = 2.;
assert(all(ismembertol(rk.r(:)', r(:)', tol, 'ByRows', true)))
assert(all(ismembertol(rk.A(:)', A(:)', tol, 'ByRows', true)))


%% Test SSP43 multistart
rk = rk_opt(4,3,'erk','ssp','startvec','smart','np',1,'write_to_file',0);
b = [1./6 1./6 1./6 1./2]';
r = 2.;
assert(all(ismembertol(rk.r(:)', r(:)', 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', 'ByRows', true)))


%% Test SSP43 multistart parallel
rk = rk_opt(4,3,'erk','ssp','startvec','smart','np',2,'write_to_file',0);
b = [1./6 1./6 1./6 1./2]';
r = 2.;
assert(all(ismembertol(rk.r(:)', r(:)', 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', 'ByRows', true)))


%% Test SDIRK SSP22 multistart
rk = rk_opt(2,2,'sdirk','ssp','startvec','random','np',1,'write_to_file',0);
b = [1./2 1./2]';
r = 4.;
assert(all(ismembertol(rk.r(:)', r(:)', 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', 'ByRows', true)))


%% Test SDIRK SSP22 multistart parallel
rk = rk_opt(2,2,'sdirk','ssp','startvec','random','np',2,'write_to_file',0);
b = [1./2 1./2]';
r = 4.;
assert(all(ismembertol(rk.r(:)', r(:)', 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', 'ByRows', true)))


%% Test DIRK SSP22 multistart
tol = 1.e-7;
rk = rk_opt(2,2,'dirk','ssp','startvec','smart','np',1,'write_to_file',0);
b = [1./2 1./2]';
r = 4.;
if abs(rk.r-r)>tol || max(max(abs(rk.b-b)))>tol
    error('Failed to find optimal DIRK SSP(2,2) method. This occasionally happens; try running the tests again.')
end
assert(all(ismembertol(rk.r(:)', r(:)', tol, 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', tol, 'ByRows', true)))


%% Test DIRK SSP22 multistart parallel
tol = 1.e-7;
rk = rk_opt(2,2,'dirk','ssp','startvec','smart','np',2,'write_to_file',0);
b = [1./2 1./2]';
r = 4.;
if abs(rk.r-r)>tol || max(max(abs(rk.b-b)))>tol
    error('Failed to find optimal DIRK SSP(2,2) method. This occasionally happens; try running the tests again.')
end
assert(all(ismembertol(rk.r(:)', r(:)', tol, 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', tol, 'ByRows', true)))


%% Test RK22 acc multistart
tol = 1.0e-6;
x = [0.117831902493812 0.187227391881097   0.812772608118903  -0.099952256513284 -0.010000000000000];
rk = rk_opt(2,2,'erk','acc','startvec',x,'np',1,'write_to_file',0);
b = [0.25 0.75]';
assert(all(ismembertol(rk.errcoeff', 1./6, 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', tol, 'ByRows', true)))


%% Test RK22 acc multistart parallel
tol = 1.0e-6;
x = [0.117831902493812 0.187227391881097   0.812772608118903  -0.099952256513284 -0.010000000000000];
rk = rk_opt(2,2,'erk','acc','startvec',x,'np',2,'write_to_file',0);
b = [0.25 0.75]';
assert(all(ismembertol(rk.errcoeff', 1./6, 'ByRows', true)))
assert(all(ismembertol(rk.b(:)', b(:)', tol, 'ByRows', true)))

%% Test RK32 3S acc
rk = rk_opt(4,2,'3Sstar','acc','algorithm','interior-point','write_to_file',0);
p = check_RK_order(rk.A,rk.b,rk.c);
g = errcoeff(rk.A,rk.b,rk.c,p);
assert(p==2)
assert(g<1.e-8)

%% Test embedded RK435 3Ss
% This is the method from table 6 of Ketcheson JCP 2010
X=[1.384996869124138 3.878155713328178 -2.324512951813145 -0.514633322274467 0.081252332929194 -1.083849060586449 -1.096110881845602 0.075152045700771 0.211361016946069 1.100713347634329 0.728537814675568 0.393172889823198 1.642598936063715 0.188295940828347 2.859440022030827 -0.655568367959557 -0.194421504490852];
s=5;
class='3Sstaremb';
[A,Ahat,b,bhat,c,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,class);
p = check_RK_order(A,b,c);
assert(p==4)
g = errcoeff(A,b,c,p);
assert(abs(g-0.005521417971326)<1.e-16);
p_hat = check_RK_order(A,bhat,c);
assert(p_hat==3)
g_hat = errcoeff(A,bhat,c,p_hat);
assert(abs(g_hat-0.063779596487680)<1.e-16);
