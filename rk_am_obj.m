function [r,g]=rk_am_obj(coeffs)
% Used by optrk

r=coeffs(end);
g=zeros(size(coeffs));
g(end)=1;
