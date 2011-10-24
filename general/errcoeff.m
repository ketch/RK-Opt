function D=errcoeff(x,class,s,p)
%function D=errcoeff(x,class,s,p)
%Computes norm of vector of leading error coefficients.
%For now we just use Butcher's approach.  We could alternatively use Albrecht's.

[A,b,c]=unpack_rk(x,s,class);


if p==1
  tau(1)=c'*b-1/2;
elseif p==2   
  tau(1)=c'.^2*b-1/3;
  tau(2)=b'*A*c-1/6;
elseif p==3
  tau(1)=c'.^3*b-1/4;
  tau(2)=(b'.*c')*A*c-1/8;
  tau(3)=b'*A*c.^2-1/12;
  tau(4)=b'*A^2*c-1/24;
elseif p==4
  tau(1)=c'.^4*b-1/5;
  tau(2)=(b.*c.^2)'*A*c-1/10;
  tau(3)=b'*(A*c).^2-1/20;
  tau(4)=(b.*c)'*A*c.^2-1/15;
  tau(5)=b'*A*c.^3-1/20;
  tau(6)=(b.*c)'*A^2*c-1/30;
  tau(7)=b'*A*diag(c)*A*c-1/40;
  tau(8)=b'*A^2*c.^2-1/60;
  tau(9)=b'*A^3*c-1/120;
elseif p==5
  tau(1)=c'.^5*b-1/6;
  tau(2)=b'*diag(c).^3*A*c-1/12;
  tau(3)=b'*diag(c)*(A*c).^2-1/24;
  tau(4)=b'*diag(c).^2*A*c.^2-1/18;
  tau(5)=b'*((A*c.^2).*(A*c))-1/36;
  tau(6)=b'*diag(c)*A*c.^3-1/24;
  tau(7)=b'*A*c.^4-1/30;
  tau(8)=b'*diag(c).^2*A^2*c-1/36;
  tau(9)=b'*((A^2*c).*(A*c))-1/72;
  tau(10)=b'*diag(c)*A*diag(c)*A*c-1/48;
  tau(11)=b'*A*diag(c).^2*A*c-1/60;
  tau(12)=b'*A*(A*c).^2-1/120;
  tau(13)=b'*diag(c)*A^2*c.^2-1/72;
  tau(14)=b'*A*diag(c)*A*c.^2-1/90;
  tau(15)=b'*A^2*c.^3-1/120;
  tau(16)=b'*diag(c)*A^3*c-1/144;
  tau(17)=b'*A*diag(c)*A^2*c-1/180;
  tau(18)=b'*A^2*diag(c)*A*c-1/240;
  tau(19)=b'*A^3*c.^2-1/360;
  tau(20)=b'*A^4*c-1/720;
elseif p==6
  disp('Order conditions for p>=6 are not coded up yet');
end

D=norm(tau);
