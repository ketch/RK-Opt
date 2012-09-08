function [r,g]=dwrk_am_obj(coeffs)
% function [r,g]=dwrk_am_obj(coeffs)
% Used by opt_dwrk

r=coeffs(end);
g=zeros(size(coeffs));
g(end)=1;
