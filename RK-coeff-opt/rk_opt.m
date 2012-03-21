function rk = rk_opt(s,p,class,objective,varargin)
%function rk = rk_opt(s,p,class,objective,varargin)
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

% Parse require, optional and param input values with inputParser
i_p = inputParser;
i_p.FunctionName = 'rk_opt';

optional_params = varargin;
setup_params(s,p,class,objective,i_p,optional_params);


% New random seed every time
rand('twister', sum(100*clock)); 

% Open pool of sessions (# is equal to the processors specified in np)
if i_p.Results.np>1
    matlabpool('local',i_p.Results.np);
end

%Set optimization parameters:
options=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-13,'MaxIter',10000,'Diagnostics','off','Display','notify','DerivativeCheck','off'...%);
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

for i=1:i_p.Results.max_tries
    fprintf('Attempt %d\n', i);

    if i_p.Results.np>1
        % # of starting points for multistart
        nsp = 20;
        n=set_n(s,class);
        
        for i=1:nsp
            x(i,:)=initial_guess(s,p,class,i_p.Results.startvec);
            tpoints = CustomStartPointSet(x);
        end

        problem = createOptimProblem('fmincon','x0',x(1,:),'objective', ...
                  @(x) rk_obj(x,class,s,p,objective),'Aeq',Aeq,'beq',beq,...
                  'lb',lb,'ub',ub,'nonlcon',...
                  @(x) nlc(x,class,s,p,objective,i_p.Results.poly_coeff_ind,...
                  i_p.Results.poly_coeff_val),'options',opts);

        ms = MultiStart('Display','final','UseParallel','always');
        [X,FVAL,status,outputg,manyminsg] = run(ms,problem,tpoints);

    else
        x=initial_guess(s,p,class,i_p.Results.startvec);

        %Optionally find a feasible (for the order conditions) point to start
        if i_p.Results.solveorderconditions==1
            x=fsolve(@(x) oc(x,class,s,p),x,opts);
        end
        
        [X,FVAL,status]=fmincon(@(x) rk_obj(x,class,s,p,objective),x,[],[],Aeq,beq,lb,ub,@(x) nlc(x,class,s,p,objective,i_p.Results.poly_coeff_ind,i_p.Results.poly_coeff_val),opts);
    end
    
    % Check order of the scheme
    if (class(1:2)=='2S' | class(1:2)=='3S')
        [rk.A,rk.b,rk.c,rk.alpha,rk.beta,rk.gamma1,rk.gamma2,rk.gamma3,rk.delta]=unpack_lsrk(X,s,class);
    else
        [rk.A,rk.b,rk.c]=unpack_rk(X,s,class);
    end
    
    order = check_RK_order(rk.A,rk.b,rk.c);
    
    % If a solution is found then exit the loop
    if (status>0 && p==order)
        break;
    end
end 

if (i==i_p.Results.max_tries && status<=0)
    fprintf('Failed to find a solution.\n')
    rk = -1;
    return
end
fprintf('The method found has order of accuracy: %d \n', order)

% Compute properties of the method
if strcmp(objective,'ssp')
    rk.r = -FVAL;
else
   % rk.r=am_radius(rk.A,rk.b,rk.c);
end
rk.errcoeff=errcoeff(rk.A,rk.b,rk.c,order);
[rk.v,rk.alpha,rk.beta] = optimal_shuosher_form(rk.A,rk.b,rk.c);
    
if (i_p.Results.writeToFile == 1 && p == order)
    output=writeFile(rk,p);
end

if i_p.Results.np>1 matlabpool close; end
end
% =========================================================================


% =========================================================================
function pop_i_p= setup_params(s,p,class,objective,i_p,optional_params)
%function pop_i_p= setup_params(s,p,class,objective,optional_params,i_p)
%
% Set default optional and param values

% Expected values or string
expected_classes = {'erk','irk','dirk','sdirk','2S','2Sstar','3Sstar'};
expected_objectives = {'ssp','acc'};
expected_startvec = {'random','smart'};
expected_solveorderconditions = [0,1];

% Default values
default_poly_coeff_ind = [];
default_poly_coeff_val = [];
default_startvec = 'random';
default_solveorderconditions = 0;
default_np = 1;
default_max_tries = 10;
default_writeToFile = 1;


% Populate input parser object
% ----------------------------
% Required values
i_p.addRequired('s',@isnumeric);
i_p.addRequired('p',@isnumeric);
i_p.addRequired('class',@(x) ischar(x) && any(validatestring(x,expected_classes)));
i_p.addRequired('objective',@(x) ischar(x) && any(validatestring(x,expected_objectives)));

% Optional values
i_p.addParamValue('poly_coeff_ind',default_poly_coeff_ind,@isnumeric); 
i_p.addParamValue('poly_coeff_val',default_poly_coeff_val,@isnumeric);
i_p.addParamValue('startvec',default_startvec,@(x) ischar(x) && any(validatestring(x,expected_startvec)));
i_p.addParamValue('solveorderconditions',default_solveorderconditions,@(x) isnumeric(x) && any(x==expected_solveorderconditions))
i_p.addParamValue('np',default_np,@isnumeric);
i_p.addParamValue('max_tries',default_max_tries,@isnumeric);
i_p.addParamValue('writeToFile',default_writeToFile,@isnumeric);

% Parse input arguments
i_p.parse(s,p,class,objective,optional_params{:});

end
% =========================================================================



% =========================================================================
function x=initial_guess(s,p,class,starttype)
%function x=initial_guess(s,p,class,starttype)
%
% Set initial guess for RK coefficients
% Includes some good initial guesses for optimal SSP methods
if ~ischar(starttype)
    x=starttype;
else
x=[];
switch class
    case 'irk'
        % Implicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s)=sort(rand(1,s)); 
                x(s+1:2*s-1)=rand(1,s-1); 
                x(2*s)=1-sum(x(s+1:2*s-1));
                x(2*s+1:2*s+s^2)=rand(1,s^2);
                x(2*s+s^2+1)=-0.01;
            case 'smart'
                x(1:s)=(1:s)/s;
                x(s+1:2*s)=1/s;
                for i=1:s
                    x(2*s+s*(i-1)+1:2*s+s*(i-1)+i-1)=1/s;
                    x(2*s+s*(i-1)+i)=1/(2*s);
                end
                x(2*s+s^2+1)=-0.01;
                x=x.*(1+rand(size(x))/4);
        end

    case 'irk5'
        % Implicit SSP methods of order >=5 always have one row of A
        % equal to zero
        switch starttype
            case 'random'
                x(1:s-1)=sort(rand(1,s-1));    %c's
                x(s:2*s-2)=rand(1,s-1);        %b's
                x(2*s-1)=1-sum(x(s:2*s-2)); 
                x(2*s:2*s-1+s*(s-1))=rand(1,s*(s-1));  %A's
                x(2*s+s*(s-1))=-0.01;            %r
            case 'smart'
                x(1:s-1)=(1:s-1)/s;
                x(s:2*s-1)=1/s;
                for i=2:s
                    x(2*s+s*(i-2):2*s+s*(i-2)+i-2)=1/s;
                    x(2*s-1+s*(i-2)+i-1)=1/(2*s);
                end
                x(2*s+s*(s-1))=-0.1;
                x=x.*(1+rand(size(x))/4);
        end

    case 'dirk'
        %Diagonally implicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s)=sort(rand(1,s));                     %c's
                x(s+1:2*s-1)=rand(1,s-1);                   %b's
                x(2*s)=1-sum(x(s+1:2*s-1));                 %last b
                x(2*s+1:2*s+s*(s+1)/2)=rand(1,s*(s+1)/2);   %A's
                x(2*s+s*(s+1)/2+1)=-0.01;                   %r
            case 'smart'
                x(1:s)=(1:s)/s;                             %c's
                x(s+1:2*s)=1/s;                             %b's
                x(2*s+1:2*s+s*(s+1)/2)=1/s;                 %A's
                j=0;
                for i=1:s
                  j=j+i;
                  x(2*s+j)=1/(2*s);                          %Diagonal A's (A_ii)
                end
                x(2*s+s*(s+1)/2+1)=-(s-(p-2)+sqrt(s^2-(p-2)))/2;                   %r
                x=x.*(1+rand(size(x))/10);
        end

    case 'sspdirk5'
        % Implicit SSP methods of order >=5 always have one row of A
        % equal to zero
        switch starttype
            case 'random'
                x(1:s-1)=sort(rand(1,s-1));    % c
                x(s:2*s-2)=rand(1,s-1)/s;      % b
                x(2*s-1)=1-sum(x(s:2*s-2));    %b(end)
                x(2*s:2*s-2+s*(s+1)/2)=rand(1,s*(s+1)/2-1)/3.;
                x(2*s+s*(s+1)/2-1)=-0.01;
            case 'smart'
                nbz=4;
                %modified 2nd order
                x(1:s-1)=(1:s-1)/s;
                x(s:2*s-1)=1/(s-nbz);  
                %Zero some b's:
                x(2*s-nbz:2*s-1)=0;
                x(s)=x(s)/2;
                x(2*s:2*s-2+s*(s+1)/2)=1/s;
                j=0;
                for i=2:s
                  j=j+i;
                  x(2*s+j-i)=1/(2*s);                          %First column A's (A_i1)
                  x(2*s+j-i+1)=3/(2*s);                        %Second column A's (A_i2)
                  x(2*s-1+j)=1/(2*s);                          %Diagonal A's (A_ii)
                end
                x(2*s+s*(s+1)/2-1)=-0.01;
                x=x.*(1+rand(size(x))/10);
        end

    case 'sdirk'
        % Singly diagonally implicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s)=sort(rand(1,s));                     %c's
                x(s+1:2*s-1)=rand(1,s-1)/s;                   %b's
                x(2*s)=1-sum(x(s+1:2*s-1));                 %last b
                x(2*s+1:2*s+1+s*(s-1)/2)=rand(1,s*(s-1)/2+1)/s;   %A's
                x(2*s+1)=rand/3/s;
                x(2*s+1+s*(s-1)/2+1)=-0.01;                   %r
        end

    case 'erk'
        %Explicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s-1)=sort(rand(1,s-1)-1/2); 
                x(s:2*s-2)=rand(1,s-1)-1/2; 
                x(2*s-1)=1-sum(x(s:2*s-2));
                x(2*s:2*s-1+s*(s-1)/2)=rand(1,s*(s-1)/2)-1/2;
                x(2*s+s*(s-1)/2)=-0.01;
            case 'smart'
                r=s-(p-3)-sqrt(s-(p-3));
                x(1:s-1)=(1:s-1)/r;
                x(s:2*s-2)=1/s;
                x(2*s-1)=1-sum(x(s:2*s-2));
                x(2*s:2*s-1+s*(s-1)/2)=1/r;
                x(2*s+s*(s-1)/2)=-r;
        end

    otherwise
        % Low-storage methods
        n=set_n(s,class);
        x=rand(1,n);
end

end
end
% =========================================================================


% =========================================================================
function wf=writeFile(rk,p)
%function wf=writeFile(rk,p)
%
% 
% Write to file Butcher's coefficients and low-storage coefficients if 
% required.

szA = size(rk.A);

outputFileName = strcat('ERK-',num2str(p),'-',num2str(szA(1)),'.txt');
writeFid = fopen(outputFileName,'w');

fprintf(writeFid, '%s\t\t %s\n', '#stage','order');
output = [szA(1);p];
fprintf(writeFid, '%u\t \t\t%u\n\n',output);

values = struct2cell(rk);
names  = fieldnames(rk);
for i=1:length(values)
    writeField(writeFid,names{i},values{i});
end

str = '==============================================================';
fprintf(writeFid,'\n%s\r\n\n',str);

wf= 1;
end
% =========================================================================
