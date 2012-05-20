function D=errcoeff(A,b,c,p)
%function D=errcoeff(A,b,c,p)
%
% Inputs: A,b,c -- Butcher tableau
%         p     -- order of accuracy of the method
% Computes the norm of the vector of truncation error coefficients
% for the terms of order p+1: 
% (elementary weight - 1/(density of the tree)/(symmetry of the tree)
% 
%
% For now we just use Butcher's approach.  We could alternatively use Albrecht's.

if p==1
  % order 2 conditions:
  tau(1) = (b'*(c) - 1/2)/1;
elseif p==2
  % order 3 conditions:
  tau(1) = (b'*((A*c)) - 1/6)/1;
  tau(2) = (b'*(c.^2) - 1/3)/2;
elseif p==3
  % order 4 conditions:
  tau(1)= (b'*((A*c.^2)) - 1/12)/2;
  tau(2)= (b'*((A*(A*c))) - 1/24)/1;
  tau(3)= (b'*((A*c).*c) - 1/8)/1;
  tau(4)= (b'*(c.^3) - 1/4)/6;
elseif p==4
  % order 5 conditions:
  tau(1) = (b'*diag(c)^3*c - 1/5)/24;
  tau(2) = (b'*diag(c)^2*A*c - 1/10)/2;
  tau(3) = (b'*diag(c)*A*diag(c)*c - 1/15)/2
  tau(4) = (b'*diag(c)*A^2*c - 1/30)/1;
  tau(5) = (b'*diag(A*c)*A*c - 1/20)/2;
  tau(6) = (b'*A*diag(c)^2*c - 1/20)/6
  tau(7) = (b'*A*diag(c)*A*c - 1/40)/1;
  tau(8) = (b'*A^2*diag(c)*c - 1/60)/2;
  tau(9) = (b'*A^3*c - 1/120)/1;
elseif p==5
  % order 6 conditions:
  tau(1) = (b'*diag(c)^4*c - 1/6)/120;
  tau(2) = (b'*diag(c)^3*A*c - 1/12)/6; 
  tau(3) = (b'*diag(c)*diag(A*c)*A*c - 1/24)/2; 
  tau(4) = (b'*diag(c)^2*A*diag(c)*c - 1/18)/4; 
  tau(5) = (b'*diag(A*c)*A*diag(c)*c - 1/36)/2; 
  tau(6) = (b'*diag(c)*A*diag(c)^2*c - 1/24)/6; 
  tau(7) = (b'*A*diag(c)^3*c - 1/30)/24;  
  tau(8) = (b'*diag(c)^2*A^2*c - 1/36)/2; 
  tau(9) = (b'*diag(A*c)*A^2*c - 1/72)/1; 
  tau(10) = (b'*diag(c)*A*diag(c)*A*c - 1/48)/1; 
  tau(11) = (b'*A*diag(c)^2*A*c - 1/60)/2; 
  tau(12) = (b'*A*diag(A*c)*A*c - 1/120)/2; 
  tau(13) = (b'*diag(c)*A^2*diag(c)*c - 1/72)/2; 
  tau(14) = (b'*A*diag(c)*A*diag(c)*c - 1/90)/2; 
  tau(15) = (b'*A^2*diag(c)^2*c - 1/120)/6; 
  tau(16) = (b'*diag(c)*A^3*c - 1/144)/1; 
  tau(17) = (b'*A*diag(c)*A^2*c - 1/180)/1; 
  tau(18) = (b'*A^2*diag(c)*A*c - 1/240)/1; 
  tau(19) = (b'*A^3*diag(c)*c - 1/360)/2; 
  tau(20) = (b'*A^4*c - 1/720)/1; 
else disp('Calculation of the principal error norm for p > 5 is not implemented yet')
end

D=norm(tau);
