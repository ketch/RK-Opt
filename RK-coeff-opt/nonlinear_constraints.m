function [con,coneq]=nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k)
% function [con,coneq]=nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k)
% Impose nonlinear constraints:
%   - if objective = 'ssp' : both order conditions and absolute monotonicity conditions
%   - if objective = 'acc' : order conditions
% The input arguments are:
%     * :math:`x`: vector of the decision variables.  See unpack_rk.m for details about
%       the order in which they are stored.
%     * *class*: class of method to search ('erk' = explicit RK; 'irk' = implicit RK; 'dirk' = diagonally implicit RK; 'sdirk' = singly diagonally implicit RK; '2S', '3S', '2S*', '3S*' = low-storage formulations).
%     * :math:`s`:number of stages.
%     * :math:`p`: order of the RK scheme.
%     * *objective*: objective function ('ssp' = maximize SSP coefficient; 'acc' = minimize leading truncation error coefficient).
%     * *poly_coeff_ind*: index of the polynomial coefficients (:math:`\beta_j`) for :math:`j > p`.
%     * *poly_coeff_val*: values of the polynomial coefficients (:math:`\beta_j`) for :math:`j > p` (tall-tree elementary weights).
% 
% The outputs are:
%     * *con*: inequality constraints, i.e. absolute monotonicity conditions if objective = 'ssp' or nothing if objective = 'acc'
%     * *coneq*: order conditions plus stability function coefficients constraints (tall-tree elementary weights)
% 
% Two forms of the order conditions are implemented: one based on **Butcher's
% approach**, and one based on **Albrecht's approach**. One or the other may lead 
% to a more tractable optimization problem in some cases, but this has not been 
% explored carefully. The Albrecht order conditions are implemented up to order 9, assuming
% a certain stage order, while the Butcher order conditions are implemented up to order 9 but
% do not assume anything about the stage order. Albrecht's approach is used
% by default.


oc_form = 'albrecht';

if k==1
    [A,b,c]=unpack_rk(x,s,class);
else
    [A,Ahat,b,bhat,D,theta] =  unpack_msrk(x,s,k,class);
end

if strcmp(objective,'ssp')
    z=-x(end); %Radius of absolute monotonicity

    if k==1 % RK methods
        %=====================================================
        % Inequality constraints: absolute monotonicity conditions
        es=ones(s,1);
        K=[A es*0;b' 0];
        G=eye(s+1)+z*K;

        con1=G\K;
        con2=G\[es;1];

        con=-[con1(:);con2(:)];
        %=====================================================

    else % multistep RK methods
        if strcmp(class(end),'2')       % Type 2 methods
            A= [zeros(k-1,s+k-1);Ahat,A];
            b= [bhat, b];
            D= [eye(k);D];
            s=length(A);
        end

        %=====================================================
        % Absolute monotonicity conditions (from Spijker 2007)
        % Construct Spijker form 
        S=zeros(s+k+1,k);
        S(1:k,1:k)=eye(k); 
        S(k+1:k+s,1:k)=D;
        S(end,1:end)=theta; 

        T=zeros(s+k+1,k+s+1);
        T(k+1:k+s,k+1:k+s)=A;
        T(k+s+1,k+1:k+s)=b;

        con1=(eye(size(T,1))+z*T)\[S z*T]; 
        con=-con1(:);
    end
else % Other objectives are handled in the objective function call
    con=[];
end

%=====================================================
% Order conditions
if k>1
    coneq= oc_ksrk(A,b,D,theta,p);
elseif strcmp(oc_form,'albrecht')
    coneq = oc_albrecht(A,b,c,p);
elseif strcmp(oc_form,'butcher')
    coneq = oc_butcher(A,b,c,p);
end
%=====================================================

for i=1:length(poly_coeff_ind)
    %Enforce stability function coefficient constraints
    j = poly_coeff_ind(i);
    coneq(end+1) = b'*A^(j-2)*c - poly_coeff_val(i);
end
