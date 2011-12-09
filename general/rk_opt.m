function rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions)
%function rk_opt(s,p,objective,restart)
% Find optimal RK methods 
%
% Variable meanings:
% s    - # of stages
% p    - order of accuracy
%
% Decision variables:
% A, b, c -- Method coefficients
% A     is s x s
% b     is s x 1
%
%
% Stored in a single vector x as:
% x=[A b' c']
% A is stored row-by-row
%
%objective:
% Set to 'ssp' to maximize SSP coefficient 
% Set to 'acc' to minimize leading truncation error coefficients
%
% poly_coeff_ind:    indices of stability polynomial coefficients to constrain
% poly_coeff_val: values of constrained stability polynomial coefficients
%
% Problem definition:
% Class of methods to search
% Available classes:
%       'erk'   : Explicit Runge-Kutta methods
%       'irk'   : Implicit Runge-Kutta methods
%       'dirk'  : Diagonally implicit Runge-Kutta methods
%       'sdirk' : Singly diagonally implicit Runge-Kutta methods
%       '2S', etc. : Low-storage explicit methods
%
%solveorderconditions:
%       if set to 1, solve the order conditions first before trying to optimize
%       (in rare cases, this is helpful for high order methods)

if nargin<8 
    solveorderconditions=0; 
end
if nargin<7 
    startvec='random';
end
if nargin<5 
    poly_coeff_ind = [];
    poly_coeff_val = [];
end

%Set optimization parameters:
options=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-13,'MaxIter',10000,'Diagnostics','off','Display','off','DerivativeCheck','off'...%);
,'Algorithm','sqp');
%For difficult cases, it can be useful to limit the line search step size
%by appending to the line above (possibly with a modified value of RelLineSrchBnd):
%,'RelLineSrchBnd',0.1,'RelLineSrchBndDuration',100000000);
%Also, sometimes something can be gained by adjusting 'Tol*' above.
%==============================================
 
if strcmp(objective,'ssp')
    opts = optimset(options,'GradObj','on');
elseif strcmp(objective,'acc')
    opts = optimset(options,'GradObj','off');
else
    error('Unrecognized objective type.');
end
                 
%Set the linear constraints: Aeq*x = beq
%and the upper and lower bounds on the unknowns: lb <= x <= ub
[Aeq,beq,lb,ub] = linear_constraints(s,class);

info=-2;
while (info==-2 || info==0) % Don't stop until we converge to a solution

if ~ischar(startvec)
    x=startvec;
else
    rand('twister', sum(100*clock)); %New random seed every time
    x=initial_guess(s,p,class,startvec);
end

%Optionally find a feasible (for the order conditions) point to start
if solveorderconditions==1
    x=fsolve(@(x) oc(x,class,s,p),x,opts);
end

% Solve the optimization problem
[X,FVAL,info]=fmincon(@(x) rk_obj(x,class,s,p,objective),x,[],[],Aeq,beq,lb,ub,@(x) nlc(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val),opts);
if strcmp(objective,'ssp')
    rk.r=-FVAL;
elseif strcmp(objective,'acc')
    rk.errcoeff=FVAL;
end

end %while loop
%==============================================

%Now extract the Butcher array of the solution from x
if class(1:2)=='2S'
    [A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,s,class)
else
    [rk.A,rk.b,rk.c]=unpack_rk(X,s,class);
end
