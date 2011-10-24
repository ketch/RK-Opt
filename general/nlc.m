function [con,coneq]=nlc(x,class,s,p)
% Nonlinear constraints for SSP RK Methods
% Including both order conditions and absolute monotonicity conditions

global objective talltree_numbers talltree_values oc_form

[A,b,c]=unpack_rk(x,s,class);

if strcmp(objective,'ssp')
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
else
    con=[];
end

%=====================================================
% Order conditions
if strcmp(oc_form,'albrecht')
    coneq = oc_albrecht(A,b,c,p);
elseif strcp(oc_form,'butcher')
    coneq = oc_butcher(A,b,c,p);
end
%=====================================================

for i=1:length(talltree_numbers)
    %Enforce stability function coefficient constraints
    j = talltree_numbers(i);
    coneq(end+1) = b'*A^(j-2)*c - talltree_values(i);
end
