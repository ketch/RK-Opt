function r=emb_obj(x,rk,phat)
%function r=emb_obj(x,rk)

%Get Butcher array
A=rk.A; b=x; c=rk.c; %bhat=x;

eps=0.;%000001;

if phat==4
  tau(1)=(c'.^4*b-1/5)/24;
  tau(2)=((b.*c.^2)'*A*c-1/10)/2;
  tau(3)=(b'*(A*c).^2-1/20)/2;
  tau(4)=((b.*c)'*A*c.^2-1/15)/2;
  tau(5)=(b'*A*c.^3-1/20)/6;
  tau(6)=((b.*c)'*A^2*c-1/30);
  tau(7)=(b'*A*diag(c)*A*c-1/40);
  tau(8)=(b'*A^2*c.^2-1/60)/2;
  tau(9)=(b'*A^3*c-1/120);
elseif phat==5
    tau(1)=(c'.^5*b-1/6)/120.;
    tau(2)=(b'*diag(c).^3*A*c-1/12)/6.;
    tau(3)=(b'*diag(c)*(A*c).^2-1/24)/2.;
    tau(4)=(b'*diag(c).^2*A*c.^2-1/18)/4.;
    tau(5)=(b'*((A*c.^2).*(A*c))-1/36)/2.;
    tau(6)=(b'*diag(c)*A*c.^3-1/24)/6.;
    tau(7)=(b'*A*c.^4-1/30)/24.;
    tau(8)=(b'*diag(c).^2*A^2*c-1/36)/2.;
    tau(9)=(b'*((A^2*c).*(A*c))-1/72)/1.;
    tau(10)=(b'*diag(c)*A*diag(c)*A*c-1/48)/1.;
    tau(11)=(b'*A*diag(c).^2*A*c-1/60)/2.;
    tau(12)=(b'*A*(A*c).^2-1/120)/2.;
    tau(13)=(b'*diag(c)*A^2*c.^2-1/72)/2.;
    tau(14)=(b'*A*diag(c)*A*c.^2-1/90)/2.;
    tau(15)=(b'*A^2*c.^3-1/120)/6.;
    tau(16)=(b'*diag(c)*A^3*c-1/144)/1.;
    tau(17)=(b'*A*diag(c)*A^2*c-1/180)/1.;
    tau(18)=(b'*A^2*diag(c)*A*c-1/240)/1.;
    tau(19)=(b'*A^3*c.^2-1/360)/2.;
    tau(20)=(b'*A^4*c-1/720)/1.;
end

%    norm(tau2)
%  tau-eps
  r= norm(abs(abs(tau)-eps))+1.e-9/min(abs(b-x));
 
 r=norm(tau);

