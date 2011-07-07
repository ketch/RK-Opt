function r=lowstorage_obj(x,class,regs)
%function [r]=lowstorage_obj(x)

%Get Butcher array
[A,b,bhat,c]=x2Abc(regs,class,x);

%r=-norm(b-bhat);
    coneq2(18)=(c'.^5*b-1/6)/120.;
    coneq2(19)=(b'*diag(c).^3*A*c-1/12)/6.;
    coneq2(20)=(b'*diag(c)*(A*c).^2-1/24)/2.;
    coneq2(21)=(b'*diag(c).^2*A*c.^2-1/18)/4.;
    coneq2(22)=(b'*((A*c.^2).*(A*c))-1/36)/2.;
    coneq2(23)=(b'*diag(c)*A*c.^3-1/24)/6.;
    coneq2(24)=(b'*A*c.^4-1/30)/24.;
    coneq2(25)=(b'*diag(c).^2*A^2*c-1/36)/2.;
    coneq2(26)=(b'*((A^2*c).*(A*c))-1/72)/1.;
    coneq2(27)=(b'*diag(c)*A*diag(c)*A*c-1/48)/1.;
    coneq2(28)=(b'*A*diag(c).^2*A*c-1/60)/2.;
    coneq2(29)=(b'*A*(A*c).^2-1/120)/2.;
    coneq2(30)=(b'*diag(c)*A^2*c.^2-1/72)/2.;
    coneq2(31)=(b'*A*diag(c)*A*c.^2-1/90)/2.;
    coneq2(32)=(b'*A^2*c.^3-1/120)/6.;
    coneq2(33)=(b'*diag(c)*A^3*c-1/144)/1.;
    coneq2(34)=(b'*A*diag(c)*A^2*c-1/180)/1.;
    coneq2(35)=(b'*A^2*diag(c)*A*c-1/240)/1.;
    coneq2(36)=(b'*A^3*c.^2-1/360)/2.;
    coneq2(37)=(b'*A^4*c-1/720)/1.;
 
 r=norm(coneq2);
