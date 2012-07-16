function [b,c] = scaled_binomial_basis(N,r,z)
% Given an arbitrary radius r, scaled_binomial generates a basis of polynomials (1+z/r)^j. 
% N is the order of polynomial basis desired.  
%
% Returns:
%  1. A matrix b, whose jth row contains the coefficients of the jth basis function
%     with the coefficients in order of ascending degree
%
%  2. A matrix c, whose jth column contains the values of the jth basis function
%     evaluated at the points z.
%
% Some of the loops could be vectorized, but since this routine isn't a bottleneck we've
% opted for clarity.

b = zeros(N+1,N+1);

for j=0:N
    for k=0:j
        b(j+1,k+1) = nchoosek(j,k)/r^k;
    end
end

if nargin < 3
    return
end

c = zeros(length(z),N+1);

for j=0:N
    for i=1:length(z)
        c(i,j+1) = (1+z(i)/r)^j;
    end
end
