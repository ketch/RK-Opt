function [X,r,errcoeff]=loop_fmincon(startvec,solveorderconditions,class,s,p,opts,objective,Aeq,beq,lb,ub,poly_coeff_ind,poly_coeff_val)
%function X=loop_fmincon(startvec,solveorderconditions,class,s,p,opts,objective,Aeq,beq,lb,ub,poly_coeff_ind,poly_coeff_val)
%
% In this function fmincon is called inside a "while loop" with some 
% suitable initial guess. The loop does not stop until it converges to a 
% solution.

info=-2; 
% info = -2: No feasible point found.
% info = 0 : Too many function evaluations or iterations.

while (info==-2 || info==0) % Don't stop until we converge to a solution

if ~ischar(startvec)
    x=startvec;
else
    rand('twister', sum(100*clock)); %New random seed every time
    x=initial_guess(s,p,class,startvec);
end

%Optionally find a feasible (for the order conditions) point to start
if solveorderconditions==1
    x=fsolve(@(x) oc(x,class,s,p),x,opts);
end

% Solve the optimization problem
[X,FVAL,info]=fmincon(@(x) rk_obj(x,class,s,p,objective),x,[],[],Aeq,beq,lb,ub,@(x) nlc(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val),opts);

end %while loop
%==============================================

% Set the objective values
if strcmp(objective,'ssp')
    r=-FVAL;
    errcoeff=[];
elseif strcmp(objective,'acc')
    errcoeff=FVAL;
    r=[];
end