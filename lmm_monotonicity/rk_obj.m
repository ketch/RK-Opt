function [r,g]=rk_obj(x,class,s,p,objective)
% function [r,g]=rk_obj(x,class,s,p,objective)
% Objective function for RK optimization.
%
% The meaning of the input arguments is as follow:
%     * :math:`x`: vector of the unknowns.
%     * class: class of method to search ('erk' = explicit RK; 'irk' = implicit RK; 'dirk' = diagonally implicit RK; 'sdirk' = singly diagonally implicit RK; '2S', '3S', '2S*', '3S*' = low-storage formulations).
%     * :math:`s`:number of stages.
%     * :math:`p`: order of the RK scheme.
%     * objective: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient).
%
% The meaning of the output arguments is as follow:
%     * r: it is a scalar containing the radius of absolute monotonicity if objective = 'ssp' or the value of the leading truncation error coefficient if objective = 'acc'.
%     * g: it is a vector and contains the gradient of the objective function respect to the unknowns.  It is an array with all zero elements except for the last component which is equal to one if objective = 'ssp' or it is an empty array if objective = 'acc'. 


if strcmp(objective,'ssp')
    r=x(end);
    g=zeros(size(x));
    g(end)=1;
elseif strcmp(objective,'acc') % Will fail for multistep-RK methods
    [A,b,c]=unpack_rk(x,s,class);
    r=errcoeff(A,b,c,p);
    g=[];
end
