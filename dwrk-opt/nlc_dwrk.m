function [con,coneq]=nlc_dwrk(x,class,s,p)
% Nonlinear constraints for SSP ARK Methods
% Including both order conditions and absolute monotonicity conditions
% Used by opt_dwrk.m

es=ones(s,1);

%=====================================================
% Extract Butcher arrays A,b,c from x
% x = [c' b' A ct' bt' At z] (A is stored row-by-row)
% Different things are stored depending on which class of methods is considered
[A,b,c,At,bt,ct]=unpack_dwrk(x,s,class);
z=-x(end); %Radius of absolute monotonicity
%=====================================================
%=====================================================
% Inequality constraints:
% Absolute monotonicity conditions
% See, e.g., Ketcheson thesis p. 34
zero_s=zeros(s,1);
K=[A  zero_s;b'  0];
Kt=[At zero_s;bt' 0];
G=eye(s+1)+z*K+z*Kt;

con1=G\K;
con2=G\Kt;
con3=G\[es;1];
con=-[con1(:);con2(:);con3(:)];
%=====================================================

%=====================================================
% Order conditions

% Underlying scheme:
A1=A-At; b1=b-bt; c1=c-ct;

% In case of 1st order -- we've already imposed this
% via a linear constraint
coneq(1)=0;

% 2nd order
if p>=2
  coneq(1)=c1'*b1-1/2.;
end
%coneq(2)=A(1,2);

if p>=3   % 3rd order
  coneq(2) =c1' .*c1' *b1 -1/3;
  coneq(3) =b1' *A1 *c1 -1/6;
end

if p>=4   % 4th order
  coneq(4)=c1'.^3*b1-1/4;
  coneq(5)=(b1'.*c1')*A1^2*es-1/8;
  coneq(6)=b1'*A1*c1.^2-1/12;
  coneq(7)=b1'*A1^2*c-1/24;
end

%HACK
%Get optimal Shu-Osher arrays as well
%r=z;
%zero_s=zeros(s,1);
%K=[A  zero_s;b'  0];
%Kt=[At zero_s;bt' 0];
%G=(eye(s+1)+r*(K+Kt));
%P=r*(G\K); Pt=r*(G\Kt); d=G\ones(s+1,1);
%Pt(3,2)=0.; coneq=[coneq(:);Pt(:)];
%END HACK

if p>=5   % 5th order
  coneq(8) =c1'.^4*b1-1/5;
  coneq(9) =(b1.*c1.^2)'*A1*c1-1/10;
  coneq(10)=b1'*(A1*c1).^2-1/20;
  coneq(11)=(b1.*c1)'*A1*c1.^2-1/15;
  coneq(12)=b1'*A1*c1.^3-1/20;
  coneq(13)=(b1.*c1)'*A1^2*c1-1/30;
  coneq(14)=b1'*A1*diag(c1)*A1*c1-1/40;
  coneq(15)=b1'*A1^2*c1.^2-1/60;
  coneq(16)=b1'*A1^3*c1-1/120;
end

if p>=6   % 6th order
  coneq(17)=c1'.^5*b1-1/6;
  coneq(18)=b1'*diag(c1).^3*A1*c1-1/12;
  coneq(19)=b1'*diag(c1)*(A1*c1).^2-1/24;
  coneq(20)=b1'*diag(c1).^2*A1*c1.^2-1/18;
  coneq(21)=b1'*((A1*c1.^2).*(A1*c1))-1/36;
  coneq(22)=b1'*diag(c1)*A1*c1.^3-1/24;
  coneq(23)=b1'*A1*c1.^4-1/30;
  coneq(24)=b1'*diag(c1).^2*A1^2*c1-1/36;
  coneq(25)=b1'*((A1^2*c1).*(A1*c1))-1/72;
  coneq(26)=b1'*diag(c1)*A1*diag(c1)*A1*c1-1/48;
  coneq(27)=b1'*A1*diag(c1).^2*A1*c1-1/60;
  coneq(28)=b1'*A1*(A1*c1).^2-1/120;
  coneq(29)=b1'*diag(c1)*A1^2*c1.^2-1/72;
  coneq(30)=b1'*A1*diag(c1)*A1*c1.^2-1/90;
  coneq(31)=b1'*A1^2*c1.^3-1/120;
  coneq(32)=b1'*diag(c1)*A1^3*c1-1/144;
  coneq(33)=b1'*A1*diag(c1)*A1^2*c1-1/180;
  coneq(34)=b1'*A1^2*diag(c1)*A1*c1-1/240;
  coneq(35)=b1'*A1^3*c1.^2-1/360;
  coneq(36)=b1'*A1^4*c1-1/720;
end

%conzero_At = At(:,1:end-1);
%conzero_At3 = At(2:3,end);
%conzero_bt = bt(1:end-1);
%conzero_A  = A( :,end);
%conzero_A2 = A(:,2)-A(1,2);
%conzero_At2 = At(:,end)-At(1,end);
%conzero_bt2 = bt(end)-At(1,end);
%conzero_b = b(end);
%%coneq = [coneq conzero_At(:)' conzero_bt' conzero_At3(:)'];
%coneq = [coneq conzero_At(:)' conzero_bt' conzero_A(:)' conzero_At2(:)' conzero_bt2 conzero_b];
