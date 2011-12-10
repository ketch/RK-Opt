function [X,r,errcoeff]=multistart_fmincon(startvec,solveorderconditions,class,s,p,opts,objective,Aeq,beq,lb,ub,poly_coeff_ind,poly_coeff_val,nsp,parallel)
%function [X,r,errcoeff]=multistart_fmincon(startvec,solveorderconditions,class,s,p,opts,objective,Aeq,beq,lb,ub,poly_coeff_ind,poly_coeff_val,nsp,parallel)
%
% In this function fmincon is called in combination with the multistart
% solver (Global Optimization Toolbox). It starts a local solver (in 
% Optimization Toolbox) from multiple starting points and stores local and 
% global solutions found during the search process.

% # of decision variables
n=set_n(s,class);

%==============================================
if ~ischar(startvec)
    for i=1:nsp
        x(i,:)=startvec;
        tpoints = CustomStartPointSet(x);
    end
else
    for i=1:nsp
        rand('twister', sum(100*clock)); %New random seed every time
        x(i,:)=initial_guess(s,p,class,startvec);
        tpoints = CustomStartPointSet(x);
    end
end
%Set initial guess
%for i=1:nsp
%    x(i,:)=rand(1,n);
%    tpoints = CustomStartPointSet(x);
%end
%==============================================

if parallel == 1
    matlabpool open
end

problem = createOptimProblem('fmincon','x0',x(1,:),'objective', ...
          @(x) rk_obj(x,class,s,p,objective),'Aeq',Aeq,'beq',beq,...
          'lb',lb,'ub',ub,'nonlcon',...
          @(x) nlc(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val), ...
          'options',opts);
ms = MultiStart('Display','final','UseParallel','always');
[X,FVAL,flagg,outputg,manyminsg] = run(ms,problem,tpoints);

if strcmp(objective,'ssp')
    r=-FVAL;
    errcoeff=[];
elseif strcmp(objective,'acc')
    errcoeff=FVAL;
    r=[];
end

if parallel == 1
    matlabpool close
end