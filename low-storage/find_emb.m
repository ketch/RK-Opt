% Find an embedded error estimator for a given method
%
% Need to put optimization in here!
%
% Variable meanings:
% s    - # of stages
% p    - order of accuracy
% regs - # of registers to allow
% class- plain, star, emb, or staremb
%
% Decision variables:
% Low-storage coefficients (fill in details later)
%

%Assumes the following are already set correctly:
% A,b,s,p

s=length(b);
rk.A=A; rk.b=b; rk.c=c; rk.m=s;
p=rk_order(rk);

maxerrcoeff=2.e-3;
rand('twister', sum(100*clock)); %New random seed every time

%Set optimization parameters:
opts=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-13,'GradObj','off','MaxIter',500,'Diagnostics','on','Display','iter');%,...
%'Algorithm','interior-point');
%==============================================

%Now set up the number of unknowns
n=s;
lb=-20*ones(n,1);
ub= 20*ones(n,1);
%==============================================

%x=rand(n,1);
%X=fsolve(@(x) nlc_emb(A,c,x,p-1),x,opts)
FVAL=10.;
while FVAL>maxerrcoeff
    %==============================================
    %Set initial guess
    x=2.*(rand(n,1)-0.5);
    %==============================================

    [X,FVAL,info]=fmincon(@(x) emb_obj(x,rk,p-1),x,[],[],[],[],lb,ub,@(x) nlc_emb(A,c,x,p-1),opts)
    %==============================================
end

%Get Butcher array
bhat=X;
rk2.A=A;rk2.b=bhat;rk2.c=c;rk2.m=s;
plotstabreg(rk2);
