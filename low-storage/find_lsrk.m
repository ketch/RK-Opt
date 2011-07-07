% Find Low-storage RK methods 
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

rand('twister', sum(100*clock)); %New random seed every time

if restart==0
    %==============================================
    %Editable options:
    regs=2; class='plain';
    s=4; p=4;
    starttype='random'; %Options: random
    maxerrcoeff=1.45e-2;
    %==============================================
end

%Set optimization parameters:
opts=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-13,'TolX',1.e-13,'GradObj','off','MaxIter',10000,'Diagnostics','on','Display','iter'...%);
,'Algorithm','interior-point');
%==============================================

%Now set up the number of unknowns
n=set_n(s,class,regs);
lb=-20*ones(n,1);
ub= 20*ones(n,1);
%==============================================

FVAL=10.;
while FVAL>maxerrcoeff
    %==============================================
    %Set initial guess
    if restart==1
      x=X;
    else
      x=rand(n,1);
    end
    %==============================================

    [X,FVAL,info]=fmincon(@(x) lowstorage_obj(x,class,regs),x,[],[],[],[],lb,ub,@(x) oc_lowstorage(regs,class,x,p),opts)
    %X=fsolve(@(x) oc_lowstorage(regs,class,x,p),x,opts)
    %==============================================
end

%Get Butcher array
[A,b,bhat,c]=x2Abc(regs,class,X);
rk.A=A;rk.b=b;rk.c=c;rk.m=s;
plotstabreg(rk);
