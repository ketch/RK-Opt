function [con,coneq]=nlc_emb(A,c,bhat,phat)
%function [con,coneq]=nlc_emb(A,c,bhat,phat)
%Order conditions

clear coneq
con=[];
%==================================================
%Order conditions for embedded method

  coneq2(1)=sum(bhat)-1;
  if phat>=2 coneq2(2)=c'*bhat-1/2; end

  if phat>=3
    coneq2(3)=c'.^2*bhat-1/3;
    coneq2(4)=bhat'*A*c-1/6;
  end

  if phat>=4
    coneq2(5)=(bhat'.*c')*A*c-1/8;
    coneq2(6)=bhat'*A*c.^2-1/12;
    coneq2(7)=bhat'*A^2*c-1/24;
    coneq2(8)=c'.^3*bhat-1/4;
  end

  if phat>=5   % 5th order
    coneq2(9) =c'.^4*bhat-1/5;
    coneq2(10) =(bhat.*c.^2)'*A*c-1/10;
    coneq2(11)=bhat'*(A*c).^2-1/20;
    coneq2(12)=(bhat.*c)'*A*c.^2-1/15;
    coneq2(13)=bhat'*A*c.^3-1/20;
    coneq2(14)=(bhat.*c)'*A^2*c-1/30;
    coneq2(15)=bhat'*A*diag(c)*A*c-1/40;
    coneq2(16)=bhat'*A^2*c.^2-1/60;
    coneq2(17)=bhat'*A^3*c-1/120;
  end

  if phat>=6   % 6th order
    coneq2(18)=c'.^5*bhat-1/6;
    coneq2(19)=bhat'*diag(c).^3*A*c-1/12;
    coneq2(20)=bhat'*diag(c)*(A*c).^2-1/24;
    coneq2(21)=bhat'*diag(c).^2*A*c.^2-1/18;
    coneq2(22)=bhat'*((A*c.^2).*(A*c))-1/36;
    coneq2(23)=bhat'*diag(c)*A*c.^3-1/24;
    coneq2(24)=bhat'*A*c.^4-1/30;
    coneq2(25)=bhat'*diag(c).^2*A^2*c-1/36;
    coneq2(26)=bhat'*((A^2*c).*(A*c))-1/72;
    coneq2(27)=bhat'*diag(c)*A*diag(c)*A*c-1/48;
    coneq2(28)=bhat'*A*diag(c).^2*A*c-1/60;
    coneq2(29)=bhat'*A*(A*c).^2-1/120;
    coneq2(30)=bhat'*diag(c)*A^2*c.^2-1/72;
    coneq2(31)=bhat'*A*diag(c)*A*c.^2-1/90;
    coneq2(32)=bhat'*A^2*c.^3-1/120;
    coneq2(33)=bhat'*diag(c)*A^3*c-1/144;
    coneq2(34)=bhat'*A*diag(c)*A^2*c-1/180;
    coneq2(35)=bhat'*A^2*diag(c)*A*c-1/240;
    coneq2(36)=bhat'*A^3*c.^2-1/360;
    coneq2(37)=bhat'*A^4*c-1/720;
  end

coneq=coneq2;
con=[];

%rk.A=A;rk.b=b;rk.c=c;rk.m=length(b);
%plotstabreg(rk);
%pause(0.01)
