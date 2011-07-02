function [con,coneq]=nlc(x,class,s,p)
% Nonlinear constraints for SSP RK Methods
% Including both order conditions and absolute monotonicity conditions
% Used by optrk

%=====================================================
es=ones(s,1);
[A,b,c]=unpack_rk(x,s,class);

z=-x(end); %Radius of absolute monotonicity
%=====================================================

%=====================================================
% Inequality constraints: absolute monotonicity conditions
K=[A es*0;b' 0];
G=eye(s+1)+z*K;

con1=G\K;
con2=G\[es;1];

con=-[con1(:);con2(:)];
%=====================================================

%=====================================================
%Equality constraints: order conditions
coneq=oc(x,class,s,p);
%=====================================================
