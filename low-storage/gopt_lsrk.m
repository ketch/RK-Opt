% Find Low-storage RK methods with small principal error norm
% Global search version
%
%Before running, type 'matlabpool open'
%
%
% Variable meanings:
% s    - # of stages
% p    - order of accuracy
% regs - # of registers to allow
% class- plain, star, emb, or staremb
%
% Decision variables:
% Low-storage coefficients (fill in details later)

clear all;
restart=0;
rand('twister', sum(100*clock)); %New random seed every time

%==============================================
%Editable options:
regs=2
class='plain';
s=8; p=5;
starttype='random'; %Options: random
maxerrcoeff=8.e-4;
nsp = 8; % # of start points for MultiSearch
%==============================================

%Set optimization parameters:
opts=optimset('MaxFunEvals',1000000,'TolCon',1.e-13,'TolFun',1.e-12,'TolX',1.e-12,'GradObj','off','MaxIter',5000,'Diagnostics','on','Display','iter'...%);
,'Algorithm','interior-point');
%==============================================

%Now set up the number of unknowns
n=set_n(s,class,regs);
lb=-20*ones(n,1);
ub= 20*ones(n,1);
%==============================================

r=10.;
while r>maxerrcoeff
    %==============================================
    %Set initial guess
    for i=1:nsp
      x(i,:)=rand(1,n);
      tpoints = CustomStartPointSet(x);
    end
    %==============================================

    problem = createOptimProblem('fmincon','x0',x(1,:),'objective',@(x) lowstorage_obj(x',class,regs),'lb',lb,'ub',ub,'nonlcon',@(x) oc_lowstorage(regs,class,x',p),'options',opts);
    ms = MultiStart('Display','final','UseParallel','always');
    [X,r,flagg,outputg,manyminsg] = run(ms,problem,tpoints);
end

X=X';
%Get Butcher array
[A,b,bhat,c]=x2Abc(regs,class,X);
rk.A=A;rk.b=b;rk.c=c;rk.m=s;
plotstabreg(rk);
