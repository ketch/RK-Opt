function [con,coneq]=nlc(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val)
% Impose nonlinear constraints:
%   - if objective = 'ssp' : both order conditions and absolute monotonicity 
%                            conditions
%   - if objective = 'acc' : order conditions

clear coneq
clear con

oc_form = 'albrecht';

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

for i=1:length(poly_coeff_ind)
    %Enforce stability function coefficient constraints
    j = poly_coeff_ind(i);
    coneq(end+1) = b'*A^(j-2)*c - poly_coeff_val(i);
end
