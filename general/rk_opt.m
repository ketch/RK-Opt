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
class='erk';


%==============================================
%Algorithmic options:
starttype='random'; 

%Set optimization parameters:
options=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-13,'MaxIter',10000,'Diagnostics','on','Display','iter','DerivativeCheck','off'...%);
,'Algorithm','sqp');

%For difficult cases, it can be useful to limit the line search step size
%by appending to the line above (possibly with a modified value of RelLineSrchBnd):
%,'RelLineSrchBnd',0.1,'RelLineSrchBndDuration',100000000);
%Also, sometimes something can be gained by adjusting 'Tol*' above.
%==============================================

% If set to 1, solve the order conditions first before trying to optimize
% (in rare cases, this is helpful for high order methods)
solveorderconditions=0;
     

% =========================================================================
% Load and read the file containing the stability polynomial coefficients
% =========================================================================
readFileName = '/Users/parsanm/Desktop/rk-stab-coeffs/results2ndSD-2D.txt';
readFid = fopen(readFileName,'r');


% =========================================================================
% Open file for writing RK Butcher table
% =========================================================================
writeFileName = '/Users/parsanm/Desktop/rk-coeffs/RKx2-2DSD.txt';
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
        talltree_values = d(6+s-fp+1:6+s)
    end


    %Set the number of decision variables
    n=set_n(s,class);

    %Set the linear constraints: Aeq*x = beq
    %and the upper and lower bounds on the unknowns: lb <= x <= ub
    [Aeq,beq,lb,ub] = linear_constraints(s,class);
    %==============================================

    info=-2;
    while (info==-2 || info==0) % Don't stop until we converge to a solution

      if strcmp(starttype,'restart')
        x=X;
      else
        x=initial_guess(s,class,starttype);
      end

      %==============================================
      %Optionally find a feasible (for the order conditions) point to start
      if solveorderconditions==1
        x=fsolve(@(x) oc(x,class,s,p),x,opts)
      end
      %==============================================
      
      if strcmp(objective,'ssp')
          obj_func = @(x) rk_obj_ssp(x,class,s,p);
          opts = optimset(options,'GradObj','on');
      else
          obj_func = @(x) rk_obj_acc(x,class,s,p);
          opts = optimset(options,'GradObj','off');
      end


      [X,FVAL,info]=fmincon(obj_func,x,[],[],Aeq,beq,lb,ub,@(x) nlc(x,class,s,p),opts);
      r=-FVAL

    end %while loop
    %==============================================
    

    %Now extract the Butcher array of the solution from x
    if class(1:2)=='2S'
        [A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,s,class)
    else
        [A,b,c]=unpack_rk(X,s,class);
    end
    
    fprintf(writeFid, '%s\t %s\n', '#stage','order');
    
    output = [s;p];
    
    fprintf(writeFid, '%u\t %u\n\n',output);
    
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
    
    str = '==============================================================';
    fprintf(writeFid,'\n%s\r\n\n',str);
       
end


fclose(readFid);
fclose(writeFid);

