function [b,c] = scaled_binomial_basis(N,h,X)
% Given an arbitrary radius h, scaled_binomial generates a basis of polynomials (1+z/h)^j. 
% N is the order of polynomial basis desired.  
%
% Returns:
%  1. A matrix b, whose jth row contains the coefficients of the jth basis function
%     with the coefficients in order of ascending degree
%
%  2. A matrix c, whose jth column contains the values of the jth basis function
%     evaluated at the points X.

basis = zeros(N+1,N+1);

for j=0:N
    for k=0:j
        basis(j+1,k+1) = nchoosek(j,k)/h^k;
    end
end

if nargin < 3
    return
end

evalX = zeros(length(X),N+1);

for j=0:N
    for i=1:length(X)
        evalX(i,j+1) = (1+X(i)/h)^j;
    end
end
