function rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,np,max_tries,writeToFile)
%function rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,np,max_tries,writeToFile)
%
% =========================================================================
% Find optimal RK methods using MATLAB's fmincon function.
% Approaches available:
%
% Optimization of a single RK method 
% ==================================
%
% - if np>1: run in parallel. fmincon is called in combination with the 
% multistart solver (Global Optimization Toolbox). It starts a local solver 
% (in Optimization Toolbox) from multiple starting points and stores local 
% and global solutions found during the search process. This approach can 
% be run in parallel.
%
% - max_tries: maximum number of fmincon/Multistart function calls.
%
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

% Set default values of the inputs arguments if they are not passed
% The implementation of this feature in Python will be much easier.
if nargin<11
    writeToFile=1; 
end
if nargin<10 
    max_tries=10;
end
if nargin<9 
    np=1;
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

if nargin<4 
    objective = 'ssp';
end

if nargin<3 
    class = 'erk';
end


rand('twister', sum(100*clock)); %New random seed every time

% Open pool of sessions (# is equal to the processors specified in np)
if np>1
    matlabpool('local',np);
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

for i=1:max_tries

    if np>1
        % # of starting points for multistart
        nsp = 20;
        n=set_n(s,class);
        
        for i=1:nsp
            x(i,:)=initial_guess(s,p,class,startvec);
            tpoints = CustomStartPointSet(x);
        end

        problem = createOptimProblem('fmincon','x0',x(1,:),'objective', ...
                  @(x) rk_obj(x,class,s,p,objective),'Aeq',Aeq,'beq',beq,...
                  'lb',lb,'ub',ub,'nonlcon',...
                  @(x) nlc(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val), ...
                  'options',opts);

        ms = MultiStart('Display','final','UseParallel','always');
        [X,FVAL,status,outputg,manyminsg] = run(ms,problem,tpoints);

    else
        x=initial_guess(s,p,class,startvec);

        %Optionally find a feasible (for the order conditions) point to start
        if solveorderconditions==1
            x=fsolve(@(x) oc(x,class,s,p),x,opts);
        end

        [X,FVAL,status]=fmincon(@(x) rk_obj(x,class,s,p,objective),x,[],[],Aeq,beq,lb,ub,@(x) nlc(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val),opts);
    end
    
    % Check order of the scheme
    if (class(1:2)=='2S' | class(1:2)=='3S')
        [rk.A,rk.b,rk.c,rk.alpha,rk.beta,rk.gamma1,rk.gamma2,rk.gamma3,rk.delta]=unpack_lsrk(X,s,class);
    else
        [rk.A,rk.b,rk.c]=unpack_rk(X,s,class);
    end
    
    order = check_RK_order(rk.A,rk.b,rk.c);
    
    % If fmincon converged to a solution and the RK scheme satisfies the 
    % order conditions get out of the loop
    if (status>0 && p==order)
        break;
    end
end 

% Print status at screen
info = print_info(status,p,order);

% Set the objective values
if strcmp(objective,'ssp')
    rk.r=-FVAL;
    rk.errcoeff=[];
elseif strcmp(objective,'acc')
    rk.errcoeff=FVAL;
    rk.r=[];
end
    
% Write output to file if required
if (writeToFile == 1 && p == order)
    output=writeFile(rk,p);
end

% Close pool sessions
if np>1
    matlabpool close;
end


