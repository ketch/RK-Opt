function [con,coneq]=oc_lowstorage(regs,class,x,p)
%function [con,coneq]=oc_lowstorage(regs,class,x,p)
%function [coneq]=oc_lowstorage(regs,class,x,p)
%Order conditions

clear coneq

%Get Butcher array
[A,b,bhat,c,alpha,beta]=x2Abc(regs,class,x);

%==================================================
%Order conditions for principal method
coneq(1)=sum(b)-1;
coneq(2)=c'*b-1/2;

if p>=3
  coneq(3)=c'.^2*b-1/3;
  coneq(4)=b'*A*c-1/6;
end

if p>=4
  coneq(5)=(b'.*c')*A*c-1/8;
  coneq(6)=b'*A*c.^2-1/12;
  coneq(7)=b'*A^2*c-1/24;
  coneq(8)=c'.^3*b-1/4;
end

if p>=5   % 5th order
  coneq(9) =c'.^4*b-1/5;
  coneq(10) =(b.*c.^2)'*A*c-1/10;
  coneq(11)=b'*(A*c).^2-1/20;
  coneq(12)=(b.*c)'*A*c.^2-1/15;
  coneq(13)=b'*A*c.^3-1/20;
  coneq(14)=(b.*c)'*A^2*c-1/30;
  coneq(15)=b'*A*diag(c)*A*c-1/40;
  coneq(16)=b'*A^2*c.^2-1/60;
  coneq(17)=b'*A^3*c-1/120;
end

if p>=6   % 6th order
  coneq(18)=c'.^5*b-1/6;
  coneq(19)=b'*diag(c).^3*A*c-1/12;
  coneq(20)=b'*diag(c)*(A*c).^2-1/24;
  coneq(21)=b'*diag(c).^2*A*c.^2-1/18;
  coneq(22)=b'*((A*c.^2).*(A*c))-1/36;
  coneq(23)=b'*diag(c)*A*c.^3-1/24;
  coneq(24)=b'*A*c.^4-1/30;
  coneq(25)=b'*diag(c).^2*A^2*c-1/36;
  coneq(26)=b'*((A^2*c).*(A*c))-1/72;
  coneq(27)=b'*diag(c)*A*diag(c)*A*c-1/48;
  coneq(28)=b'*A*diag(c).^2*A*c-1/60;
  coneq(29)=b'*A*(A*c).^2-1/120;
  coneq(30)=b'*diag(c)*A^2*c.^2-1/72;
  coneq(31)=b'*A*diag(c)*A*c.^2-1/90;
  coneq(32)=b'*A^2*c.^3-1/120;
  coneq(33)=b'*diag(c)*A^3*c-1/144;
  coneq(34)=b'*A*diag(c)*A^2*c-1/180;
  coneq(35)=b'*A^2*diag(c)*A*c-1/240;
  coneq(36)=b'*A^3*c.^2-1/360;
  coneq(37)=b'*A^4*c-1/720;
end
%==================================================

%==================================================
%Order conditions for embedded method
if strcmp(class,'emb') || strcmp(class,'staremb')
  phat=p-1;

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
else coneq2=[];
end
  %==================================================
  %Constrain tall trees

  %Dormand-Prince (7 stages):
%  coneq(18)=b'*A^4*c-1./600.;
%  coneq(19)=b'*A^5*c;
  %DP embedded method:
%  coneq2(9)=bhat'*A^3*c-0.009141666666666666;
%  coneq2(10)=bhat'*A^4*c-0.001341666666666666;
%  coneq2(11)=bhat'*A^5*c-1./24000.;

  %Fehlberg (6 stages):
%  coneq(18)=b'*A^5*c-1./2080.;
  %Embedded method:
%  coneq2(9)=bhat'*A^4*c-1./104.;
%  coneq2(10)=bhat'*A^5*c;

  %Kennedy's RK5(4)8[3R]C
  %coneq(18)=b'*A^4*c-1./640.;
  %coneq(19)=b'*A^5*c-1./3400.;
  %coneq(20)=b'*A^6*c-2.608022e-5;

  %DDAS47:
  %coneq(end+1)=b'*A^3*c-0.00808984135736739751;
  %coneq(end+1)=b'*A^4*c-0.00114500144821215386;
  %coneq(end+1)=b'*A^5*c-0.0000933161720160215541;

  %Calvo46:
  %coneq(end+1)=b'*A^3*c-0.00785333;
  %coneq(end+1)=b'*A^4*c-0.00094889;

  %LDD46:
  %coneq(end+1)=b'*A^3*c-0.00781005;
  %coneq(end+1)=b'*A^4*c-0.00132141;

  %BB26:
  %coneq(end+1)=b'*A  *c-0.165919771368;
  %coneq(end+1)=b'*A^2*c-0.040919732041;
  %coneq(end+1)=b'*A^3*c-0.007555704391;
  %coneq(end+1)=b'*A^4*c-0.000891421261;
  %==================================================

  %Simplifying condition C(2):
  %m=length(b);
  %coneq(end+1:end+m)=A*c-c.^2/2.;

  %Force fake stages:
%  coneq(end+1)=beta(9,8);
%  coneq(end+1)=beta(10,9);
coneq=[coneq coneq2];%sum(abs(b-bhat))-0.01];
con=[];

%rk.A=A;rk.b=b;rk.c=c;rk.m=length(b);
%plotstabreg(rk);
%pause(0.01)
