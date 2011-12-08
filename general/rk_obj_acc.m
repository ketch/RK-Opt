function [r]=rk_obj_acc(x,class,s,p)
% function [r]=rk_obj_acc(x,class,s,p)
%
% Compute norm of vector of leading error coefficients
% No gradient is computed

r = errcoeff(x,class,s,p);