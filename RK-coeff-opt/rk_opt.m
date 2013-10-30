function rk = rk_opt(s,p,class,objective,varargin)
% function rk = rk_opt(s,p,class,objective,varargin)
%
% Find optimal RK and multistep RK methods.
% The meaning of the arguments is as follows:
%
%     * `s` number of stages.
%     * `k` number of steps (1 for RK methods)
%     * `p` order of the Runge-Kutta (RK) scheme.
%     * class: class of method to search.  Available classes:
%
%       * 'erk'      : Explicit Runge-Kutta methods
%       * 'irk'      : Implicit Runge-Kutta methods
%       * 'dirk'     : Diagonally implicit Runge-Kutta methods
%       * 'sdirk'    : Singly diagonally implicit Runge-Kutta methods
%       * '2S', etc. : Low-storage explicit methods; see *Ketcheson, "Runge-Kutta methods with minimum storage implementations". J. Comput. Phys. 229(5):1763 - 1773, 2010*)
%       * 'emsrk1/2'    : Explicit multistep-Runge-Kutta methods
%       * 'imsrk1/2'    : Implicit multistep-Runge-Kutta methods
%       * 'dimsrk1/2'   : Diagonally implicit multistep-Runge-Kutta methods
%
%     * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient)
%       Accuracy optimization is not currently supported for multistep RK methods
%     * poly_coeff_ind: index of the polynomial coefficients to constrain (`\beta_j`) for `j > p`  (j denotes the index of the stage). The default value is an empty array.  Note that one should not include any indices `i \le p`, since those are determined by the order conditions.
%     * poly_coeff_val: constrained values of the polynomial coefficients (`\beta_j`) for `j > p` (tall-tree elementary weights). The default value is an empty array.
%     * startvec: vector of the initial guess ('random' = random approach; 'smart' = smart approach; alternatively, the user can provide the startvec array. By default startvec is initialize with random numbers.
%     * solveorderconditions: if set to 1, solve the order conditions first before trying to optimize. The default value is 0.
%     * np: number of processor to use. If np `> 1` the MATLAB global optimization toolbox *Multistart* is used. The default value is 1 (just one core).
%     * max_tries: maximum number of fmincon function calls. The default value is 10.
%     * writeToFile: whether to write to a file. If set to 1 write the RK coefficients to a file called "ERK-p-s.txt". The default value is 1.
%     * algorithm: which algorithm to use in fmincon: 'sqp','interior-point', or 'active-set'. By default sqp is used.
%
%     .. note::
%        **numerical experiments have shown that when the objective function is the minimization of the leading truncation error coefficient, the interior-point algorithm performs much better than the sqp one.**
%
%     * display: level of display of fmincon solver ('off', 'iter', 'notify' or 'final'). The default value is 'notify'.
%     * problem_class: class of problems for which the RK is designed ('linear' or 'nonlinear' problems). This option changes the type of order conditions check, i.e. linear or nonlinear order conditions controll. The default value is 'nonlinear'.
%
%
%     .. note::
%
%        Only `s` , `p` , class and objective are required inputs.
%        All the other arguments are **parameter name - value arguments to the input
%        parser scheme**. Therefore they can be specified in any order.
%
%    **Example**::
%
%     >> rk=rk_opt(4,3,'erk','acc','max_tries',2,'np',1,'solveorderconditions',1)
%
% The fmincon options are set through the **optimset** that creates/alters optimization options structure. By default the following additional options are used:
%     * MaxFunEvals = 1000000
%     * TolCon = 1.e-13
%     * TolFun = 1.e-13
%     * TolX = 1.e-13
%     * MaxIter = 10000
%     * Diagnostics = off
%     * DerivativeCheck = off
%     * GradObj = on, if the objective is set equal to 'ssp'


[k,np,max_tries,startvec,poly_coeff_ind,poly_coeff_val,...
    solveorderconditions,write_to_file,algorithm,display,min_amrad]=...
    setup_params(varargin);

% New random seed every time
rand('twister', sum(100*clock));

% Open pool of sessions (# is equal to the processors specified in np)
if np>1
    matlabpool('local',np);
end

%Set optimization parameters:
options=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-15,'MaxIter',10000,'Diagnostics','off','Display',display,'DerivativeCheck','off'...%);
,'Algorithm',algorithm);
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
[Aeq,beq,lb,ub] = linear_constraints(s,class,objective,k);

for i=1:max_tries
    fprintf('Attempt %d\n', i);

    if np>1
        % # of starting points for multistart
        nsp = 20;

        for i=1:nsp
            x(i,:) = initial_guess(s,p,class,startvec,k);
            tpoints = CustomStartPointSet(x);
        end

        problem = createOptimProblem('fmincon','x0',x(1,:),'objective', ...
                  @(x) rk_obj(x,class,s,p,objective),'Aeq',Aeq,'beq',beq,...
                  'lb',lb,'ub',ub,'nonlcon',...
                  @(x) nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,...
                  poly_coeff_val,k),'options',opts);

        ms = MultiStart('Display','final','UseParallel','always');
        [X,FVAL,status,outputg,manyminsg] = run(ms,problem,tpoints);

    else
        x = initial_guess(s,p,class,startvec,k);

        %Optionally find a feasible (for the order conditions) point to start
        if solveorderconditions==1
            x=fsolve(@(x) order_conditions(x,class,s,p,Aeq,beq),x);
        end

        [X,FVAL,status] = fmincon(@(x) rk_obj(x,class,s,p,objective),...
                                x,[],[],Aeq,beq,lb,ub,...
                                @(x) nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k),opts);
    end


    % Check order of the scheme
    if (class(1:2)=='2S' | class(1:2)=='3S')
        [rk.A,rk.Ahat,rk.b,rk.bhat,rk.c,rk.chat,rk.alpha,rk.beta,rk.gamma1,rk.gamma2,rk.gamma3,rk.delta] = unpack_lsrk(X,s,class);
        order = check_RK_order(rk.A,rk.b,rk.c,'nonlinear');
    elseif k==1
        [rk.A,rk.b,rk.c] = unpack_rk(X,s,class);
        order = check_RK_order(rk.A,rk.b,rk.c,'nonlinear');
    else
        [rk.A,rk.Ahat,rk.b,rk.bhat,rk.rk.D,rk.theta] = unpack_msrk(X,s,k,class);
        order = p; %HACK
    end

    % Compute properties of the method
    if strcmp(objective,'ssp')
        rk.r = -FVAL;
    else
        rk.r = am_radius(rk.A,rk.b,rk.c);
    end

     % If a solution is found, then exit the loop
    if (status>0 && p==order && (~strcmp(objective,'ssp') || rk.r>min_amrad))
        fprintf('The method found has order of accuracy: %d \n', order)
        break;
    end
end

if (i==max_tries && status<=0)
    fprintf('Failed to find a solution.\n')
    rk = -1;
    return
end

if k==1
    rk.errcoeff = errcoeff(rk.A,rk.b,rk.c,order);

    if strcmp(objective,'ssp')
        [rk.v_opt,rk.alpha_opt,rk.beta_opt] = optimal_shuosher_form(rk.A,rk.b,rk.c);
    end

    if (write_to_file == 1 && p == order)
        output = write_file(rk,p);
    end
end

if np>1 matlabpool close; end
end
% =========================================================================


% =========================================================================

function [k,np,max_tries,startvec,poly_coeff_ind,poly_coeff_val,...
    solveorderconditions,write_to_file,algorithm,display,min_amrad]=...
    setup_params(optional_params)
%function [k,np,max_tries,startvec,poly_coeff_ind,poly_coeff_val,...
%    solveorderconditions,write_to_file,algorithm,display,min_amrad]=...
%    setup_params(optional_params)
%
% Set default optional and param values

i_p = inputParser;
i_p.FunctionName = 'setup_params';

expected_solveorderconditions = [0,1];
expected_algorithms = {'sqp', 'interior-point','active-set'};
expected_displays = {'notify', 'iter', 'final'};
expected_problem_class = {'linear', 'nonlinear'};

% Default values
default_startvec = 'random';
default_solveorderconditions = 0;
default_np = 1;
default_max_tries = 10;
default_write_to_file = 1;
default_algorithm = 'sqp';
default_display = 'notify';
default_problem_class = 'nonlinear';


% Populate input parser object
% ----------------------------
% Parameter values
i_p.addParamValue('k',1,@isnumeric);
i_p.addParamValue('min_amrad',0,@isnumeric);
i_p.addParamValue('poly_coeff_ind',[],@isnumeric);
i_p.addParamValue('poly_coeff_val',[],@isnumeric);
i_p.addParamValue('startvec',default_startvec);
i_p.addParamValue('solveorderconditions',default_solveorderconditions,@(x) isnumeric(x) && any(x==expected_solveorderconditions))
i_p.addParamValue('np',default_np,@isnumeric);
i_p.addParamValue('max_tries',default_max_tries,@isnumeric);
i_p.addParamValue('write_to_file',default_write_to_file,@isnumeric);
i_p.addParamValue('algorithm',default_algorithm,@(x) ischar(x) && any(validatestring(x,expected_algorithms)));
i_p.addParamValue('display',default_display,@(x) ischar(x) && any(validatestring(x,expected_displays)));


i_p.parse(optional_params{:});

k                    = i_p.Results.k;
min_amrad            = i_p.Results.min_amrad;
np                   = i_p.Results.np;
max_tries            = i_p.Results.max_tries;
startvec             = i_p.Results.startvec;
poly_coeff_ind       = i_p.Results.poly_coeff_ind;
poly_coeff_val       = i_p.Results.poly_coeff_val;
solveorderconditions = i_p.Results.solveorderconditions;
write_to_file        = i_p.Results.write_to_file;
algorithm            = i_p.Results.algorithm;
display              = i_p.Results.display;
end
% =========================================================================



% =========================================================================
function x=initial_guess(s,p,class,starttype,k)
% function x=initial_guess(s,p,class,starttype,k)
%
% Set initial guess for RK coefficients
% Includes some good initial guesses for optimal SSP methods
if ~ischar(starttype)
    x=starttype;
else
x=[];
switch class
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

    otherwise
        % Low-storage and multistep-RK methods
        n=set_n(s,class,k);
        x=rand(1,n);
end

end
end
% =========================================================================


% =========================================================================
function wf=write_file(rk,p)
%function wf=write_file(rk,p)
%
%
% Write to file Butcher's coefficients and low-storage coefficients if
% required.

szA = size(rk.A);

output_file_name = strcat('ERK-',num2str(p),'-',num2str(szA(1)),'.txt');
write_fid = fopen(output_file_name,'w');

fprintf(write_fid, '%s\t\t %s\n', '#stage','order');
output = [szA(1);p];
fprintf(write_fid, '%u\t \t\t%u\n\n',output);

values = struct2cell(rk);
names  = fieldnames(rk);
for i=1:length(values)
    write_field(write_fid,names{i},values{i});
end

str = '==============================================================';
fprintf(write_fid,'\n%s\r\n\n',str);

wf= 1;
end
% =========================================================================
