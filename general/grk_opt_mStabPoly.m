%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find optimal Runge-Kutta coefficients by using the MULTISTART toolbox 
%
%
% Variable meanings:
% s                : # of stages
% p                : order of accuracy
% talltree_numbers : stability function coefficient constraints (positions) 
% talltree_values  : stability function coefficient constraints (values)  
%
% Decision variables:
% A, b, c -- RK method coefficients
% A     is s x s
% b     is s x 1
%
%
% Stored in a single vector x as:
% x = [A b' c']
% A is stored row-by-row
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restart optimization procedure
restart = 0;

if restart==0
    clear all;
    restart = 0;
end

% New random seed every time
rand('twister', sum(100*clock)); 


% Define gloabl variables
global oc_form objective talltree_numbers talltree_values

       

% ==============================================
% Problem definition:
% ==============================================

% Order condition form
oc_form = 'albrecht'  % Choices: albrecht or butcher

% Objective
objective = 'acc'     % Set to 'ssp' to maximize SSP coefficient 
                      % Set to 'acc' to minimize the leading truncation 
                      % error coefficients


% Class of methods to search
% Available classes:
%       'erk'      : Explicit Runge-Kutta methods
%       'irk'      : Implicit Runge-Kutta methods
%       'dirk'     : Diagonally implicit Runge-Kutta methods
%       'sdirk'    : Singly diagonally implicit Runge-Kutta methods
%       '2S', etc. : Low-storage explicit methods
class='erk';



% ==============================================
% Algorithmic options:
% ==============================================
starttype = 'random';

options = optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-16,...
                   'TolX',1.e-16,'GradObj','off','MaxIter',5000,...
                   'Diagnostics','on','Display','iter','Algorithm','sqp');
               
               
% Whether to solve order condition first before trying to optimize
% In rare cases, this is helpful for high order methods
solveorderconditions = 0;

% # of start points for MultiSearch
nsp = 5; 


% =========================================================================
% Open file for writing RK coefficients
% =========================================================================
outputFileName = strcat('results/ERK-x2.txt');
writeFid = fopen(outputFileName,'w'); 

               

% =========================================================================
% Load and read the file containing the stability polynomial coefficients
% =========================================================================
readFileName = 'ADD';
readFid = fopen(readFileName,'r');

% Header
tline = fgets(readFid); 

% Read number of stability polynomial written in the file 
tline = fgets(readFid);

% Convert string to floating numbers
d = sscanf(tline,'%f');
nbrStabPoly = sscanf(tline,'%f');
nbrStabPoly = floor(nbrStabPoly);

% Read White line
tline = fgets(readFid);
 
% Header
tline = fgets(readFid);

% Read White line
tline = fgets(readFid);


% Loop over the stability polynomial
for i_stabPoly = 1:1
    
    % Read information
    tline = fgets(readFid);
    
    % Convert string to floating numbers
    d = sscanf(tline,'%f');
    
    % Number of stages:
    s = floor(d(1));
    
    % Order of accuracy:
    p = floor(d(2));
    
    % Free parameters in the stability function optimization
    fp = floor(d(3));
    
    % Set tall tree numbers (indices) and tall tree values
    if fp == 0
        talltree_numbers = []
        talltree_values = []
    else
        talltree_numbers = s-fp+1:s
        talltree_values = d(6+s-fp+2:length(d))
    end
    
    % Define objective function
    if strcmp(objective,'ssp')
        obj_fun = @(x) rk_obj_ssp(x,class,s,p);
        opts = optimset(options,'GradObj','on');
    else
        obj_fun = @(x) rk_obj_acc(x,class,s,p);
        opts = optimset(options,'GradObj','off');
    end
    
    
    % ==============================================
    % Set up the problem and solve it:
    % ==============================================
    
    % # number of decision variables
    n = set_n(s,class);
    
    % Set the linear constraints: Aeq*x = beq
    % and the upper and lower bounds on the unknowns: lb <= x <= ub
    [Aeq,beq,lb,ub] = linear_constraints(s,class,n);
    
    % Set initial guess
    for i=1:nsp
        x(i,:) = min(lb) + (max(ub) - min(lb)).*rand(1,n); %rand(1,n);
        tpoints = CustomStartPointSet(x);
    end
    
    % Set the nonlinear constraits class
    % oc_fh is the order conditions function handle
    if strcmp(oc_form,'albrecht')
        oc_fh = @oc_albrecht;
    elseif strcmp(oc_form,'butcher')
        oc_fh = @oc_butcher;
    end
    nonlcon = @(x) nlc(x,class,s,p,oc_fh);
    
    % Define the optimization problem
    problem = createOptimProblem('fmincon','x0',x(1,:),'objective',obj_fun,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',opts);
    
    % Create multistart object
    ms = MultiStart('UseParallel','always','Display','iter');
    
    % Run the multistart optimization algorithm
    [X,fval,exitflag,output,solutions] = run(ms,problem,tpoints)
    
    
    
    % ==============================================
    % Write coefficients to file:
    % ==============================================
    
    % Open file for writing RK Butcher table
    % outputFileName = strcat('results/','ERK-',num2str(p),'-',num2str(s),'.txt');
    % writeFid = fopen(outputFileName,'w');  
    % fprintf(writeFid, '%s\t %s\n', '#stage','order');
    output = [s;p];
    fprintf(writeFid, '%u\t %u\n\n',output);
    
    % Now extract the Butcher array of the solution from x and the low storage
    % coefficients if required
    if (class(1:2)=='2S' | class(1:2)=='3S')
        [A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta] = unpack_lsrk(X,s,class)
        
        % Write to file Butcher's coefficients
        str = 'A';
        if writeFid ~= -1
            fprintf(writeFid,'%s\r\n',str);
        end
        [rows cols] = size(A);
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],A');
        
        str = 'b';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(b');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],b');
        
        str = 'c^T';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(c');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n\n'],c');
        
        % Write to file coefficients for low storage formulation
        str = 'alpha';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(alpha');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],alpha');
        
        str = 'beta';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(beta');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],beta');
        
        str = 'gamma1';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(gamma1');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],gamma1');
        
        str = 'gamma2';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(gamma2');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],gamma2');
        
        str = 'gamma3';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(gamma3');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],gamma3');
        
        str = 'delta';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(delta');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],delta');
        
        
        str = '==============================================================';
        fprintf(writeFid,'\n%s\r\n\n',str);
    else
        [A,b,c] = unpack_rk(X,s,class);
        
        % Write to file Butcher's coefficients
        str = 'A';
        if writeFid ~= -1
            fprintf(writeFid,'%s\r\n',str);
        end
        [rows cols] = size(A);
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],A');
        
        str = 'b';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(b');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n'],b');
        
        str = 'c^T';
        fprintf(writeFid,'\n%s\r\n',str);
        [rows cols] = size(c');
        x = repmat('%5.16E\t',1,(cols-1));
        fprintf(writeFid,[x,'%5.16E\n\n'],c');
        
    end
    
    clear x
    clear s
    clear talltree_numbers
    clear talltree_values
    clear n
    
end


fclose(readFid);
fclose(writeFid);



