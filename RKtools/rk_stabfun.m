function [p,q] = rk_stabfun(rk)
% function [p,q] = rk_stabfun(rk)
%
% Outputs the stability function of a RK method.
% The Butcher coefficients should be stored in rk.A, rk.b.
%
% p contains the coefficients of the numerator
%
% q contains the coefficients of the denominator
%
% $$\\phi(z)=\\frac{\\sum_j p_j z^j}{\\sum_j q_j z^j} = \\frac{\\det(I-z(A+eb^T))}{\\det(I-zA)}.$$

s = length(rk.b);
e=ones(s,1);
p=poly(rk.A-e*rk.b'); %Numerator
q=poly(rk.A);         %Denominator
