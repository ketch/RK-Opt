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
                  
%Set optimization parameters:
options=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-13,'MaxIter',10000,'Diagnostics','on','Display','iter','DerivativeCheck','off'...%);
,'Algorithm','sqp');

%For difficult cases, it can be useful to limit the line search step size
%by appending to the line above (possibly with a modified value of RelLineSrchBnd):
%,'RelLineSrchBnd',0.1,'RelLineSrchBndDuration',100000000);
%Also, sometimes something can be gained by adjusting 'Tol*' above.
%==============================================
 
if strcmp(objective,'ssp')
    obj_fun = @(x) rk_obj_ssp(x,class,s,p);
    opts = optimset(options,'GradObj','on');
else
    obj_fun = @(x) rk_obj_acc(x,class,s,p);
    opts = optimset(options,'GradObj','off');
end
                 
talltree_numbers=[2]
talltree_values=[0.5]

%==============================================
% Problem definition:
% Class of methods to search
% Available classes:
%       'erk'   : Explicit Runge-Kutta methods
%       'irk'   : Implicit Runge-Kutta methods
%       'dirk'  : Diagonally implicit Runge-Kutta methods
%       'sdirk' : Singly diagonally implicit Runge-Kutta methods
%       '2S', etc. : Low-storage explicit methods
class='2S';

%Number of stages:
s=4; 

%Order of accuracy:
p=1;

%==============================================
%Algorithmic options:
starttype='random'; 

%if set to 1, solve the order conditions first before trying to optimize
%(in rare cases, this is helpful for high order methods)
solveorderconditions=0;

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

  [X,FVAL,info]=fmincon(@(x) rk_obj(x,class,s,p),x,[],[],Aeq,beq,lb,ub,@(x) nlc(x,class,s,p),opts);
  r=-FVAL

end %while loop
%==============================================

%Now extract the Butcher array of the solution from x
if class(1:2)=='2S'
    [A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,s,class)
else
    [A,b,c]=unpack_rk(X,s,class);
end
