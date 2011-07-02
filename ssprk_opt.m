% Find optimal SSP RK methods 
%
% Variable meanings:
% s    - # of stages
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

rand('twister', sum(100*clock)); %New random seed every time

%==============================================
%Editable options:
solveorderconditions=0;
class='erk';
s=9; p=4;
starttype='random'; %Options: random
%==============================================

%Set optimization parameters:
opts=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-13,'GradObj','on','MaxIter',10000,'Diagnostics','on','Display','iter','DerivativeCheck','on'...%);
,'Algorithm','sqp');
%,'Algorithm','interior-point');
%,'RelLineSrchBnd',0.1,'RelLineSrchBndDuration',100000000);
%For difficult cases, it can be useful to limit the line search step size
%by appending to the line above (with possibly modified value of RelLineSrchBnd):
%Also, sometimes something can be gained by adjusting 'Tol*' above.
%==============================================

%Now set up:
%The number of unknowns -                     n
%The linear constraints -                     Aeq*x = beq
%The upper and lower bounds on the unknowns - ub, lb

n=set_n(s,class);

[Aeq,beq,lb,ub] = linear_constraints(s,class);
%==============================================

info=-2;
while (info==-2 || info==0)
  %==============================================
  %Set initial guess
  if restart==1
    x=X;
  else
    x=initial_guess(s,class,starttype);
  end
  %==============================================

  %==============================================
  %Optionally find a feasible (for the order conditions) point to start
  if solveorderconditions==1
    x=fsolve(@(x) oc(x,class,s,p),x,opts)
  end
  %==============================================
  %The optimization call:
  [X,FVAL,info]=fmincon(@(x) rk_am_obj(x),x,[],[],Aeq,beq,lb,ub,@(x) nlc(x,class,s,p),opts)
  r=-FVAL;
end %while loop
%==============================================

%Now extract the Butcher array of the solution from x
[A,b,c]=unpack_rk(X,s,class)
