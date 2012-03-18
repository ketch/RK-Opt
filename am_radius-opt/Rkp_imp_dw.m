function [R,alpha,beta,tbeta]=Rkp_imp_dw(k,p)
%function [R,alpha,beta]=Rkp_imp_dw(k,p)
%Author: David Ketcheson
%
%Finds the optimal k-step implicit LMM with order of accuracy p
%allowing downwinding
%
%Inputs: k = # of steps
%        p = order of accuracy
%Outputs: alpha, beta, tbeta = the coefficients of the method
    
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

M=3*k+2;              %Number of decision variables (coefficients)
rmax=M+.0001; rmin=0; %Upper and lower bounds for R
r=rmax;               %Initial guess
c=zeros(M,1); d=zeros(p+1,1); B=zeros(p+1,M);
lb=zeros(M,1); ub=zeros(M,1)+1.e6;
%=========================================================

 
while (rmax-rmin>acc) %Find R by bisection
 %Set up equality constraints
 %g: First k+1 unkowns are beta's, next k+1 are \tilde{\beta}'s, last k are gammas
  for i=0:p
    d(i+1)=k^i;
    for j=0:k-1
      if (i+j==0) %Avoid divide-by-zero
        B(i+1,  j+1)=r;                       %beta
        B(i+1,k+1+j+1)=r;                    %tbeta
      else
        B(i+1,  j+1)=r*j^i + i*j^(i-1);       %beta
        B(i+1,k+1+j+1)=r*j^i - i*j^(i-1);    %tbeta
      end %if
      B(i+1,2*(k+1)+j+1)=j^i;                 %gamma
    end
    B(i+1,k+1)=i*k^(i-1);
    B(i+1,2*(k+1))=i*k^(i-1);
  end
  %Test feasibility for this value of r
  [g,lambda,exitflag]=linprog(c,[],[],B,d,lb,ub,c,opts);
  if exitflag==1
    rmin=r; r=(r+rmax)/2;
  else
    rmax=r; r=(rmin+r)/2;
  end
end

%Now get a feasible solution so we have the coefficients
r=rmin;
for i=0:p
  d(i+1)=k^i;
  for j=0:k-1
    if (i+j==0) %Avoid divide-by-zero
      B(i+1,j+1)=r;
      B(i+1,k+1+j+1)=r;
    else
      B(i+1,j+1)=r*j^i + i*j^(i-1);
      B(i+1,k+1+j+1)=r*j^i - i*j^(i-1);
    end %if
    B(i+1,2*(k+1)+j+1)=j^i;
  end
  B(i+1,k+1)=i*k^(i-1);
  B(i+1,2*(k+1))=i*k^(i-1);
end
[g,lambda,exitflag]=linprog(c,[],[],B,d,lb,ub,c,opts);
R=r;
beta=g(1:k+1); tbeta=g(k+2:2*k+2); alpha=g(2*k+3:end)+r*(beta(1:k)+tbeta(1:k));
alpha=wrev(alpha); beta=wrev(beta); tbeta=wrev(tbeta);
