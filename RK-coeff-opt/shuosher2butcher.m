function [A,b,c]=shuosher2butcher(alpha,beta)
% function [A,b,c]=shuosher2butcher(alpha,beta);
%
% Generate Butcher form of a Runge-Kutta method,
% given its Shu-Osher or modified Shu-Osher form
%
% For an m-stage method, alpha and beta (or lambda and mu) should be 
% matrices of dimension (m+1) x m 
%
% Note that MATLAB indexes from 1, while the Shu-Osher coefficients
% are usually indexed from zero.

s=size(alpha,2);
X=eye(s)-alpha(1:end-1,:);
A=X\beta(1:end-1,:);
b=beta(end,:)+alpha(end,:)*A; b=b';
c=sum(A,2);
