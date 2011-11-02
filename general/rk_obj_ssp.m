function [r,g]=rk_obj_ssp(x,class,s,p)
%function [r,g]=rk_obj_ssp(x,class,s,p)

r = x(end);

g = zeros(size(x));

g(end) = 1;