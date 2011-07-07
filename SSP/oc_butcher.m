function coneq=oc_butcher(A,b,c,p)
% function coneq=oc_butcher(A,b,c,p)
% Order conditions for RKMs
% Assumes p>1
% This version is based on Butcher's approach.

coneq(1)=c'*b-1/2;

if p>=3   
  coneq(2)=c'.^2*b-1/3;
  coneq(3)=b'*A*c-1/6;
end

if p>=4  
  coneq(4)=c'.^3*b-1/4;
  coneq(5)=(b'.*c')*A*c-1/8;
  coneq(6)=b'*A*c.^2-1/12;
  coneq(7)=b'*A^2*c-1/24;
end

if p>=5 
  coneq(8)=c'.^4*b-1/5;
  coneq(9)=(b.*c.^2)'*A*c-1/10;
  coneq(10)=b'*(A*c).^2-1/20;
  coneq(11)=(b.*c)'*A*c.^2-1/15;
  coneq(12)=b'*A*c.^3-1/20;
  coneq(13)=(b.*c)'*A^2*c-1/30;
  coneq(14)=b'*A*diag(c)*A*c-1/40;
  coneq(15)=b'*A^2*c.^2-1/60;
  coneq(16)=b'*A^3*c-1/120;
end

if p>=6
  coneq(17)=c'.^5*b-1/6;
  coneq(18)=b'*diag(c).^3*A*c-1/12;
  coneq(19)=b'*diag(c)*(A*c).^2-1/24;
  coneq(20)=b'*diag(c).^2*A*c.^2-1/18;
  coneq(21)=b'*((A*c.^2).*(A*c))-1/36;
  coneq(22)=b'*diag(c)*A*c.^3-1/24;
  coneq(23)=b'*A*c.^4-1/30;
  coneq(24)=b'*diag(c).^2*A^2*c-1/36;
  coneq(25)=b'*((A^2*c).*(A*c))-1/72;
  coneq(26)=b'*diag(c)*A*diag(c)*A*c-1/48;
  coneq(27)=b'*A*diag(c).^2*A*c-1/60;
  coneq(28)=b'*A*(A*c).^2-1/120;
  coneq(29)=b'*diag(c)*A^2*c.^2-1/72;
  coneq(30)=b'*A*diag(c)*A*c.^2-1/90;
  coneq(31)=b'*A^2*c.^3-1/120;
  coneq(32)=b'*diag(c)*A^3*c-1/144;
  coneq(33)=b'*A*diag(c)*A^2*c-1/180;
  coneq(34)=b'*A^2*diag(c)*A*c-1/240;
  coneq(35)=b'*A^3*c.^2-1/360;
  coneq(36)=b'*A^4*c-1/720;
end

if p>=7
  disp('Order conditions for p>6 are not coded up yet');
end
