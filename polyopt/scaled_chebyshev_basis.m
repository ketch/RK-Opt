function [b,c] = scaled_chebyshev_basis(N,zmin,zmax,z)
% Given an arbitrary domain on the real axis [zmin,zmax], scaled_chebyshev generates a basis of Chebyshev Polynomials of
% the first kind scaled and shifted by the affine mapping: m(x)=m1*x+m0 where m1=2/(zmax-zmin) and m0=-(1+zmin) so m([zmin,zmax]) -> [-1,1]
% N is the order of polynomial basis desired.  
%
% Returns:
%  1. A matrix b, whose ith row contains the coefficients of the ith basis function
%     with the coefficients in order of ascending degree
%
%  2. A matrix c, whose ith column contains the values of the ith basis function
%     evaluated at the points z.
%
% Some of the loops could be vectorized, but since this routine isn't a bottleneck we've
% opted for clarity.

b = zeros(N+1,N+1);

m1 = 2/(zmax-zmin);
m0 = -(1+zmin*m1);

% T_0' = 1
b(1,1) = 1;

% T_1' = m1*x + m0
b(2,1) = m0;
b(2,2) = m1;

% T_{n+1}' = 2*(m1*x + m0)*T_{n} - T_{n-1} 
for k=1:N-1 
    b(k+2,:) = 2*(m1*[0 b(k+1,1:end-1)] + m0*b(k+1,:)) - b(k,:);
end

if nargin < 4
    return
end

c = zeros(length(z),N+1);

c(:,1) = 1;
c(:,2) = m1*z+m0;

for k=1:N-1 
    c(:,k+2) = 2*(m1*c(:,k+1).*z + m0*c(:,k+1)) - c(:,k);
end
