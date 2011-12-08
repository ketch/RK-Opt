function [con,coneq]=nlc(x,class,s,p,oc_fh)
% function [con,coneq]=nlc(x,class,s,p)
%
% Nonlinear constraints for
% - SSP RK methods (absolute monotonicity condition)
% - order conditions: Butcher or Albrecht approaches
% - stability function coefficient constraints

% Global variables
global objective talltree_numbers talltree_values oc_form

% Unpack RK coefficients
[A,b,c]=unpack_rk(x,s,class);

if strcmp(objective,'ssp')
    %=====================================================
    % Inequality constraints: absolute monotonicity conditions
    es = ones(s,1);
    z = -x(end);
    K = [A es*0;b' 0];
    G = eye(s+1)+z*K;

    con1 = G\K;
    con2=  G\[es;1];

    con = -[con1(:);con2(:)];
    %=====================================================
else
    con = [];
end

%=====================================================
% Order conditions
coneq = oc_fh(A,b,c,p);
%=====================================================



%=====================================================
% Enforce stability function coefficient constraints
for i = 1:length(talltree_numbers)
    j = talltree_numbers(i);
    coneq(end + 1) = b'*A^(j-2)*c - talltree_values(i);
end
