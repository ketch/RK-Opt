function [R,gamma]=Rskp(s,k,p)
%function [R,gamma]=Rskp(s,k,p)
%Author: David Ketcheson
%
%Finds the optimal contractive k-step, s-stage GLM with order of accuracy p
%for linear problems
%
%Inputs: s = # of stages
%        k = # of steps
%        p = order of accuracy
%Outputs: 
%        R = threshold factor
%        gamma = coefficients of the polynomials
%         
%        for k=1, the resulting polynomial is
%        \sum_{j=0}^m (1+z/R)^j
%
%        in general, the resulting stability function is
%        ...
%Depends on MATLAB's optimization toolbox for the LP solver
    
%=========================================================
%Initialize
tbtest=which('linprog');
if ~exist(tbtest)
  disp('MATLAB optimization toolbox is required')
  return
end
%Set options for linprog
opts=optimset('TolX',1.e-15,'TolFun',1.e-15,'MaxIter',1000000,...
               'LargeScale','off','Simplex','off','Display','off');
acc=1.e-15; %Accuracy of bisection search

M=k*(s+1);
rmax=s+0.0001;
rmin=0;
r=rmax;  %Initial guess
c=zeros(M,1); d=zeros(p+1,1); B=zeros(p+1,M);
%=========================================================

%=========================================================
%Find R by bisection
while (rmax-rmin>acc) 
  %Set up and improve conditioning of equality constraints
  for i=0:p
    d(i+1)=k^i;
    for m=1:M
      ii=ceil(m/(s+1));
      j=m-(ii-1)*(s+1)-1;
      B(i+1,m)=0;
      for l=0:min(i,j)
        B(i+1,m)=B(i+1,m)+factorial(i)/factorial(l)/factorial(i-l) ...
          *(k-ii)^(i-l)/r^l*prod(j-(0:(l-1)));
      end
    end
  end
  %Test feasibility for this value of r
  [gamma,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts); 
  if exitflag==1
    rmin=r; r=(r+rmax)/2;
  else
    rmax=r; r=(rmin+r)/2;
  end
end
%=========================================================

%=========================================================
%Now get a feasible solution so we have the coefficients
r=rmin;
for i=0:p
  d(i+1)=k^i;
  for m=1:M
    ii=ceil(m/(s+1));
    j=m-(ii-1)*(s+1)-1;
    B(i+1,m)=0;
    for l=0:min(i,j)
      B(i+1,m)=B(i+1,m)+factorial(i)/factorial(l)/factorial(i-l) ...
        *(k-ii)^(i-l)/r^l*prod(j-(0:(l-1)));
    end
  end
end
[gamma,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts); 
%=========================================================

R=r;
