function [R,alpha,beta]=Rkp_imp(k,p)
%function [R,alpha,beta]=Rkp_imp(k,p)
%Author: David Ketcheson
%
%Find the optimal SSP k-step implicit LMM with order of accuracy p
%
%Inputs: k = # of steps
%        p = order of accuracy
%
%Outputs: alpha, beta = the coefficients of the method
%
%Depends on MATLAB's optimization toolbox for the LP solver
    
%=========================================================
%Initialize
clear tbtest;
tbtest=which('linprog');
if ~exist(tbtest)
  disp('MATLAB optimization toolbox is required')
  return
end
%Set options for linprog
opts=optimset('TolX',1.e-15,'TolFun',1.e-15,'MaxIter',10000000,...
               'LargeScale','off','Simplex','off','Display','off');
acc=1.e-15; %Accuracy of bisection search

M=2*k+1;              %Number of decision variables (coefficients)
rmax=2.0001; rmin=0;  %Upper and lower bounds for R
r=rmax;               %Initial guess
c=zeros(M,1); d=zeros(p+1,1); B=zeros(p+1,M);
%=========================================================
  
%=========================================================
%Find R by bisection
while (rmax-rmin>acc) 
  %Set up equality constraints
  %g: First k+1 unkowns are beta's, last k are gammas
  for i=0:p
    d(i+1)=1;
    for j=0:k-1
      if (i+j==0) %Avoid divide-by-zero
        B(i+1,j+1)=r/k^i;
      else
        B(i+1,j+1)=(r*j + i)*(j/k)^(i-1) / k;
      end %if
      B(i+1,k+j+2)=(j/k)^i;
    end
    B(i+1,k+1)=i/k;
  end
  %Test feasibility for this value of r
  [x,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts);
  if exitflag==1
    rmin=r; r=(r+rmax)/2;
  else
    rmax=r; r=(rmin+r)/2;
  end
end

%=========================================================
%Now get a feasible solution so we have the coefficients
r=rmin;
for i=0:p
  d(i+1)=1;
  for j=0:k-1
    if (i+j==0) %Avoid divide-by-zero
      B(i+1,j+1)=r/k^i;
    else
      B(i+1,j+1)=(r*j + i)*(j/k)^(i-1) / k;
    end %if
    B(i+1,k+j+2)=(j/k)^i;
  end
  B(i+1,k+1)=i/k;
end
[g,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts);
%=========================================================

%Prepare outputs
R=r;
beta=g(1:k+1); alpha=g(k+2:end)+r*beta(1:end-1);
