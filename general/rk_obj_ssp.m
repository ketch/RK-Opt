function [r,g]=rk_obj_ssp(x,class,s,p)
% function [r,g]= rk_obj_ssp(x,class,s,p)
% 
% Set the objective function for the optimization of the radius of
% monotonicity of the SSP scheme.
%
% x: decision variables
%    The first n-1 values are the RK coefficients
%
% r: radius of monotonicity (= x(end))
%
% g: gradient (Jacobian matrix)

r = x(end);
g = zeros(size(x));
g(end) = 1;
