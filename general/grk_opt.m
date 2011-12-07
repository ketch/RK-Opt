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

% Restart
restart=0;

if restart==0
    clear all;
    restart=0;
end

global oc_form objective talltree_numbers talltree_values

rand('twister', sum(100*clock)); %New random seed every time

oc_form = 'albrecht'  % Choices: albrecht or butcher
objective = 'acc' % Set to 'ssp' to maximize SSP coefficient 
                  % Set to 'acc' to minimize leading truncation error coefficients
                  
%==============================================
% Problem definition:
% Class of methods to search
% Available classes:
%       'erk'   : Explicit Runge-Kutta methods
%       'irk'   : Implicit Runge-Kutta methods
%       'dirk'  : Diagonally implicit Runge-Kutta methods
%       'sdirk' : Singly diagonally implicit Runge-Kutta methods
%       '2S', etc. : Low-storage explicit methods
class='3Sstar';


%==============================================
%Algorithmic options:
starttype='random'; 

%Set optimization parameters:
options=optimset('MaxFunEvals',1000000,'TolCon',1.e-10,'TolFun',1.e-10,'TolX',1.e-10,'MaxIter',10000,'Diagnostics','on','Display','iter','DerivativeCheck','off'...%);
,'Algorithm','sqp');

%For difficult cases, it can be useful to limit the line search step size
%by appending to the line above (possibly with a modified value of RelLineSrchBnd):
%,'RelLineSrchBnd',0.1,'RelLineSrchBndDuration',100000000);
%Also, sometimes something can be gained by adjusting 'Tol*' above.
%==============================================

% If set to 1, solve the order conditions first before trying to optimize
% (in rare cases, this is helpful for high order methods)
solveorderconditions=0;

% Number of starting points
nsp = 40;
     

% =========================================================================
% Load and read the file containing the stability polynomial coefficients
% =========================================================================
readFileName = '/Users/parsanm/Desktop/SD-2.txt';
readFid = fopen(readFileName,'r')



% =========================================================================
% Open file for writing RK Butcher table
% =========================================================================
writeFileName = '/Users/parsanm/Desktop/optRK-coefs/RKx2-2DSD.txt';
writeFid = fopen(writeFileName,'w');

% Header
tline=fgets(readFid); 

% Read number of stability polynomial written in the file 
tline=fgets(readFid);

% Convert string to floating numbers
d = sscanf(tline,'%f');
nbrStabPoly = sscanf(tline,'%f');
nbrStabPoly = floor(nbrStabPoly);

% Read White line
tline=fgets(readFid);

% Header
tline=fgets(readFid);

% Read White line
tline=fgets(readFid);

for i_stabPoly = 1:1
    
    % Read information
    tline=fgets(readFid);

    % Convert string to floating numbers
    d = sscanf(tline,'%f');

    % Number of stages:
    s = floor(d(1));

    % Order of accuracy:
    p = floor(d(2));

    % Free parameters in the stability function optimization
    fp = floor(d(3));

    if fp == 0
        talltree_numbers = []
        talltree_values = []
    else
        talltree_numbers = s-fp+1:s
        talltree_values = d(6+s-fp+2:length(d))
    end


    %Set the number of decision variables
    n=set_n(s,class);

    %Set the linear constraints: Aeq*x = beq
    %and the upper and lower bounds on the unknowns: lb <= x <= ub
    [Aeq,beq,lb,ub] = linear_constraints(s,class);
    %==============================================
    
    
    %maxerrcoeff=8.e-4;
    %r=10.;
    %while r>maxerrcoeff % Don't stop until we converge to a solution

    
    %==============================================
    %Set initial guess
    for i=1:nsp
        x(i,:)= rand(1,n);
        tpoints = CustomStartPointSet(x);
    end
    
    %==============================================
    
    %==============================================
    %Optionally find a feasible (for the order conditions) point to start
    if solveorderconditions==1
        x=fsolve(@(x) oc(x,class,s,p),x)
    end
    %==============================================
    
    if strcmp(objective,'ssp')
        obj_func = @(x) rk_obj_ssp(x,class,s,p);
        opts = optimset(options,'GradObj','on');
    else
        obj_func = @(x) rk_obj_acc(x,class,s,p);
        opts = optimset(options,'GradObj','off');
    end
    
    
    
    problem = createOptimProblem('fmincon','x0',x(1,:),'objective',obj_func,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',@(x) nlc(x,class,s,p),'options',opts);
    ms = MultiStart('UseParallel','always','Display','iter');
    
    %matlabpool open 8

    [X,r,flagg,outputg,manyminsg] = run(ms,problem,tpoints)
    
    %matlabpool close


    %end %while loop
    %==============================================
        

    %Now extract the Butcher array of the solution from x
    if (class(1:2)=='2S' | class(1:2)=='3S')
        [A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,s,class)
    else
        [A,b,c]=unpack_rk(X,s,class);
    end
    
    fprintf(writeFid, '%s\t %s\n', '#stage','order');
    
    output = [s;p];
    
    fprintf(writeFid, '%u\t %u\n\n',output);
    
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
    
    
    clear x
    clear s
    clear talltree_numbers
    clear talltree_values
    clear n
       
end


fclose(readFid);
fclose(writeFid);