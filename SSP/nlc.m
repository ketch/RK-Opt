function [con,coneq]=nlc(x,class,s,p)
% Nonlinear constraints for SSP RK Methods
% Including both order conditions and absolute monotonicity conditions

[A,b,c]=unpack_rk(x,s,class);

%=====================================================
% Inequality constraints: absolute monotonicity conditions
es=ones(s,1);
z=-x(end);
K=[A es*0;b' 0];
G=eye(s+1)+z*K;

con1=G\K;
con2=G\[es;1];

con=-[con1(:);con2(:)];
%=====================================================

%=====================================================
%Equality constraints: order conditions
% Would be nice to test effectiveness of Butcher and Albrecht formulations
% against each other.
%coneq=oc_albrecht(A,b,c,p);
coneq=oc_butcher(A,b,c,p);
%=====================================================
