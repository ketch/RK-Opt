%Find optimal SSP DARK methods with s stages and order p
%DARK = Downwind Additive Runge-Kutta
%Uses rational formulation, as documented in Ketcheson et. al. 2007
%Variable meanings:
% s         - # of stages
% p         - order of accuracy
% class     - type of ARK method (explicit, implicit, etc.)
%       Right now only explicit is working
% starttype - type of initial guess (random, smart, load)
%
% Decision variables:
% A,At,b,bt,c,ct - coefficients of the ARK method
%
% At, bt, ct are short for A-tilde, etc.
% The method is:
% y_i = y_{n-1} + h sum_j a_{ij} f(y_j) + sum_j at_{ij} ft(y_j)
% Where ft() is the negative downwind version of f()
%
% r     - the radius of absolute monotonicity (multiplied by -1)
% x     - the decision variables stored in the following order:
% x=[c' b' A ct' bt' At r] % (A,At are stored row-by-row)

if restart==0
  clear all;  restart=0;
  rand('twister', sum(100*clock)); %New random seed every time
end

%==============================================
%Editable options:
class='irk'; s=3; p=3;
starttype='random'; 
coeff_bound=10.;
%==============================================
tol=1.e-13;
opts=optimset('MaxFunEvals',1000000,'TolCon',tol,'TolFun',tol,'TolX',tol,'GradObj','on','MaxIter',50,'Diagnostics','on','Display','iter','Algorithm','sqp');
%,'RelLineSrchBnd',0.00000001,'RelLineSrchBndDuration',100000000);
%For difficult cases, it can be useful to limit the line search step size
%by appending to the line above (with possibly modified value of RelLineSrchBnd):
%Also, sometimes something can be gained by adjusting 'Tol*' above.

%==============================================
%Now set up:
%The number of unknowns -                     n
%The linear constraints -                     Aeq*x = beq
%The upper and lower bounds on the unknowns - ub, lb
[n,n2]=set_n_dwrk(s,class);
switch class
  case 'erk' % Explicit (A, At are strictly lower triangular)
    Aeq=zeros(2*s-1,n);
    beq=zeros(2*s-1,1);
    for i=2:s         
      Aeq(i-1,i-1)=1;     %Require c's to be row sums of A
      Aeq(i-1,2*s+(i-2)*(i-1)/2:2*s-1+i*(i-1)/2)=-1;
      Aeq(i-2+s,i-1+n2)=1;     %Require ct's to be row sums of At
      Aeq(i-2+s,2*s+(i-2)*(i-1)/2+n2:2*s-1+i*(i-1)/2+n2)=-1;
    end
    % sum(b-bt)=1
    Aeq(2*s-1,s:2*s-1)=1; Aeq(2*s-1,(s+n2):(2*s-1+n2))=-1; beq(2*s-1)=1;

  case 'dirk' % Diagonally implicit (A, At are lower triangular)
    Aeq=zeros(2*s+1,n);
    beq=zeros(2*s+1,1);
    for i=1:s   %Require c's to be row sums of A and ct's to be row sums of At      
      Aeq(i,i)=1;
      Aeq(i,(2*s+i*(i+1)/2-(i-1)):(2*s+i*(i+1)/2)) = -1;
      Aeq(i+s,i+n2)=1;
      Aeq(i+s,(2*s+i*(i+1)/2-(i-1)+n2):(2*s+i*(i+1)/2+n2)) = -1;
    end
    % sum(b-bt)=1; same code here as for irk case below
    Aeq(2*s+1,s+1:2*s)=1;  % b
    Aeq(2*s+1,(s+1+n2):(2*s+n2)) = -1;  % bt
    beq(end)=1;

  case 'irk'  %Fully implicit
    Aeq=zeros(2*s+1,n); beq=zeros(2*s+1,1);
    for i=1:s  %Require c's to be row sums of A and ct's to be row sums of At
      Aeq(i,i)=1;
      Aeq(i,((2+(i-1))*s+1):((2+i)*s))=-1;
      Aeq(i+s,i+n2)=1;
      Aeq(i+s,((2+(i-1))*s+1+n2):((2+i)*s)+n2)=-1;
    end
    % sum(b-bt)=1
    Aeq(2*s+1,s+1:2*s)=1; Aeq(2*s+1,n2+s+1:n2+2*s)=-1; beq(end)=1;  
    %ub(s+1:2*s)=1;
end
lb=-coeff_bound+zeros(1,n); lb(end)=-300;
ub=coeff_bound+zeros(1,n); ub(end)=0;
%==============================================

info=-2; FVAL=0;
%while (info==-2 || info==0)
while FVAL>-20

  if restart==1
    x=X;
  else
    x=initial_guess_dwrk(s,class,starttype);
  end

  %==============================================
  %The optimization call:
  [X,FVAL,info]=fmincon(@dwrk_am_obj,x,[],[],Aeq,beq,lb,ub,@(x) nlc_dwrk(x,class,s,p),opts)
  r=-FVAL; 
  %==============================================
  [A,b,c,At,bt,ct]=unpack_dwrk(X,s,class);
end
%==============================================

%Get optimal Shu-Osher arrays as well
r
r=-X(end)
zero_s=zeros(s,1);
K=[A  zero_s;b'  0];
Kt=[At zero_s;bt' 0];
G=(eye(s+1)+r*(K+Kt));
P=r*(G\K); Pt=r*(G\Kt); d=G\ones(s+1,1);
%==============================================


%  case 'dirk' %Diagonally Implicit (A is lower triangular)
%    n=2*s+s*(s+1)/2+1;
%    Aeq=zeros(s+1,n);
%    beq=zeros(1,s+1);
%    for i=1:s
%      Aeq(i,i)=1;
%      Aeq(i,2*s+1+i*(i-1)/2:2*s+i*(i+1)/2)=-1;
%    end
%    Aeq(end,s+1:2*s)=1; beq(end)=1;
%    lb=zeros(1,n); lb(end)=-30;
%    ub=2+zeros(1,n); ub(end)=0; ub(s+1:2*s)=1;
%
%  case 'sdirk' %Singly Diagonally Implicit
%    n=2*s+s*(s+1)/2+1-s+1;
%    Aeq=zeros(s+1,n);
%    beq=zeros(1,s+1);
%    Aeq(1,1)=1; Aeq(1,2*s+1)=-1;
%    for i=2:s
%      Aeq(i,i)=1;
%      Aeq(i,2*s+1)=-1;
%      Aeq(i,2*s+2+(i-1)*(i-2)/2:2*s+1+i*(i-1)/2)=-1;
%    end
%    Aeq(end,s+1:2*s)=1; beq(end)=1;
%    lb=zeros(1,n); lb(end)=-30;
%    ub=2+zeros(1,n); ub(end)=0; ub(s+1:2*s)=1;

%  case 'dirk'
%    switch starttype
%      case 'random'
%        x(1:s)=sort(rand(1,s));                     %c's
%        x(s+1:2*s-1)=rand(1,s-1);                   %b's
%        x(2*s)=1-sum(x(s+1:2*s-1));                 %last b
%        x(2*s+1:2*s+s*(s+1)/2)=rand(1,s*(s+1)/2);   %A's
%        x(2*s+s*(s+1)/2+1)=-0.01;                   %r
%      case 'smart'
%        x(1:s)=(1:s)/s;                             %c's
%        x(s+1:2*s)=1/s;                             %b's
%        x(2*s+1:2*s+s*(s+1)/2)=1/s;                 %A's
%	j=0;
%	for i=1:s
%          j=j+i;
%          x(2*s+j)=1/(2*s);                          %Diagonal A's (A_ii)
%	end
%        x(2*s+s*(s+1)/2+1)=-(s-(p-2)+sqrt(s^2-(p-2)))/2;                   %r
%	x=x.*(1+rand(size(x))/10);
%      case 'load'
%	x=loadX(A,b,c,r,'dirk');
%    end
%
%  case 'sdirk'
%    switch starttype
%      case 'random'
%        x(1:s)=sort(rand(1,s));                     %c's
%        x(s+1:2*s-1)=rand(1,s-1)/s;                   %b's
%        x(2*s)=1-sum(x(s+1:2*s-1));                 %last b
%        x(2*s+1:2*s+1+s*(s-1)/2)=rand(1,s*(s-1)/2+1)/s;   %A's
%	x(2*s+1)=rand/3/s;
%        x(2*s+1+s*(s-1)/2+1)=-0.01;                   %r
%  end


%    case 'sdirk'
%      c=X(1:s); b=X(s+1:2*s);
%      A(1,1)=X(2*s+1);
%      for i=2:s
%        A(i,1:i-1)=X(2*s+2+(i-2)*(i-1)/2:2*s+1+i*(i-1)/2);
%        A(i,i)=X(2*s+1);
%      end
