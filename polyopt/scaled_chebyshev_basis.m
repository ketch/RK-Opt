function [basis,evalX] = scaled_chebyshev_basis(N,a,b,X)
% Given an arbitrary domain on the real axis [a,b], scaled_chebyshev generates a basis of Chebyshev Polynomials of
% the first kind scaled and shifted by the affine mapping: m(x)=m1*x+m0 where m1=2/(b-a) and m0=-(1+a) so m([a,b]) -> [-1,1]
% N is the order of polynomial basis desired.  
%
% Returns:
%  1. A matrix basis, whose ith row contains the coefficients of the ith basis function
%     with the coefficients in order of ascending degree
%
%  2. A matrix evalX, whose ith column contains the values of the ith basis function
%     evaluated at the points X.

basis = zeros(N+1,N+1);

m1 = 2/(b-a);
m0 = -(1+a*m1);

% T_0' = 1
basis(1,1) = 1;

% T_1' = m1*x + m0
basis(2,1) = m0;
basis(2,2) = m1;

% T_{n+1}' = 2*(m1*x + m0)*T_{n} - T_{n-1} 
for k=1:N-1 
    basis(k+2,:) = 2*(m1*[0 basis(k+1,1:end-1)] + m0*basis(k+1,:)) - basis(k,:);
end

if nargin < 4
    return
end

evalX = zeros(length(X),N+1);

evalX(:,1) = 1;
evalX(:,2) = m1*X+m0;

for k=1:N-1 
    evalX(:,k+2) = 2*(m1*evalX(:,k+1).*X + m0*evalX(:,k+1)) - evalX(:,k);
end
