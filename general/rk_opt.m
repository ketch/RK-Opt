function rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,multi_start,parallel,writeToFile)
%function rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,multi_start,parallel,writeToFile)
%
% =========================================================================
% Find optimal RK methods using MATLAB's fmincon function.
% Approaches available:
%
% Optimization of a single RK method 
% ==================================
%
% - whie loop: fmincon is called inside a "while loop" with some suitable
% initial guess. The loop does not stop until it converges to a solution.
%
% - multistart solver: fmincon is called in combination with the multistart
% solver (Global Optimization Toolbox). It starts a local solver (in 
% Optimization Toolbox) from multiple starting points and stores local and 
% global solutions found during the search process. This approach can be 
% run in parallel.
% =========================================================================
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
% poly_coeff_ind: indices of stability polynomial coefficients to constrain
% poly_coeff_val: values of constrained stability polynomial coefficients
%
% Problem definition:
% Class of methods to search
% Available classes:
%       'erk'      : Explicit Runge-Kutta methods
%       'irk'      : Implicit Runge-Kutta methods
%       'dirk'     : Diagonally implicit Runge-Kutta methods
%       'sdirk'    : Singly diagonally implicit Runge-Kutta methods
%       '2S', etc. : Low-storage explicit methods
%
%solveorderconditions:
%       if set to 1, solve the order conditions first before trying to optimize
%       (in rare cases, this is helpful for high order methods)

if nargin<11
    writeToFile=0; 
end
if nargin<10 
    parallel=0; 
end
if nargin<9 
    multi_start=0; 
end
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

if multi_start==1
    % # of starting points for multistart
    nsp = 40;
    [X,rk.r,rk.errcoeff] = multistart_fmincon(startvec,solveorderconditions,class,s,p,opts,objective,Aeq,beq,lb,ub,poly_coeff_ind,poly_coeff_val,nsp,parallel); 
elseif multi_start==0
    [X,rk.r,rk.errcoeff]=loop_fmincon(startvec,solveorderconditions,class,s,p,opts,objective,Aeq,beq,lb,ub,poly_coeff_ind,poly_coeff_val);
else
    error('multi_start must be 0 or 1');
end

    
%Now extract the Butcher array of the solution from sol and the low-storage
%coefficients if request by the class.
%Then write them to a file.
if (class(1:2)=='2S' | class(1:2)=='3S')
    [rk.A,rk.b,rk.c,rk.alpha,rk.beta,rk.gamma1,rk.gamma2,rk.gamma3,rk.delta]=unpack_lsrk(X,s,class);
    if writeToFile == 1
        output=writeFile(rk.A,rk.b,rk.c,rk.alpha,rk.beta,rk.gamma1,rk.gamma2,rk.gamma3,rk.delta,s,p);
    end
else
    [rk.A,rk.b,rk.c]=unpack_rk(X,s,class);
    if writeToFile == 1
        output=writeFile(rk.A,rk.b,rk.c,[],[],[],[],[],[],s,p);
    end
end


% Simple check: order of the scheme
order = check_RK_order(rk.A,rk.b,rk.c);
if order == p
    fprintf('The RK coefficients satisfy the order conditions. \n')
    fprintf('Order of accuracy: %d \n\n', order)
else
    fprintf('The RK coefficients does not satisfy the order conditions. \n');
    fprintf('Order of accuracy: %d \n\n', order)
end


