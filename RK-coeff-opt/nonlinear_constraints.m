function [con,coneq]=nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k,emb_poly_coeff_ind,emb_poly_coeff_val,constrain_emb_stability)
% function [con,coneq]=nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k,emb_poly_coeff_ind,emb_poly_coeff_val,constrain_emb_stability)
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
%     * :math:`k`: Number of steps for multi-step, mlti-stage schemes.
%     * *emb_poly_coeff_ind*: index of the polynomial coefficients of the embedded scheme (:math:`\beta_j`) for :math:`j > p`.
%     * *emb_poly_coeff_val*: values of the polynomial coefficients of the embedded scheme (:math:`\beta_j`) for :math:`j > p` (tall-tree elementary weights).
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
    [A,b,c,Ahat,bhat,chat] = unpack_rk(x,s,class);
else
    [A,Ahat,b,bhat,D,theta] = unpack_msrk(x,s,k,class);
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
    if ~isempty(bhat)
        coneq2 = oc_albrecht(Ahat,bhat,chat,p-1);
        coneq = [coneq  coneq2];
    end
elseif strcmp(oc_form,'butcher')
    coneq = oc_butcher(A,b,c,p);
    if ~isempty(bhat)
        coneq2 = oc_butcher(Ahat,bhat,chat,p-1);
        coneq = [coneq  coneq2];
    end
end
%=====================================================

for i=1:length(poly_coeff_ind)
    %Enforce stability function coefficient constraints
    j = poly_coeff_ind(i);
    coneq(end+1) = b'*A^(j-2)*c - poly_coeff_val(i);
end

for i=1:length(emb_poly_coeff_ind)
    %Enforce stability function coefficient constraints
    j = emb_poly_coeff_ind(i);
    coneq(end+1) = bhat'*Ahat^(j-2)*chat - emb_poly_coeff_val(i);
end
%=====================================================

if ~isempty(constrain_emb_stability)
    rk_tmp.A = Ahat; rk_tmp.b = bhat; rk_tmp.c = chat;
    % matlab stores polynomial coefficients for polyval etc. in another order
    poly_coef = rk_stabfun(rk_tmp);
    poly_coef = poly_coef(end:-1:1);
    res = polyval(poly_coef, constrain_emb_stability);
    con = [con, res .* conj(res) - 1];
end

% Energy Preserving conditions
% Choose pseudo-energy-preserving q (up to 6)
q=6;
if q>=1
    coneq(end+1) = sum(b)-1;
end
if q>=2
    coneq(end+1) = c'*b-1/2;
end
if q >= 3
     coneq(end+1) = c'.^2*b-1/3;
end
if q >= 4
    coneq(end+1) = b'*A^2*c-b'*A*c+1/8;
    coneq(end+1) = c'.^3*b-1/4;
    coneq(end+1) = ((b'.*c')*A*c)-(b'*A*c.^2)/2-1/12;
end
if q >= 5
    coneq(end+1) = (b'*A^2*c.^2)/2 + (b.*c)'*A^2*c-b'*A^2*c-(b'*A*c.^2)/2 + (b'*A*c)/2 - 1/24;
    coneq(end+1) = c'.^4*b-1/5;
    coneq(end+1) = 2*(b'*A*diag(c)*A*c)-(b'*(A*c).^2) - (b'*A^2*c)-(b'*A*c.^2)+(b'*A*c)-1/24;
    coneq(end+1) = (b.*c.^2)'*A*c-(b'*A*c.^3)/3 - 1/12;
end

if q >= 6
    coneq(end+1) = -(b'*A*c)^2 + (b'*A^2*c)- 2 *(b'*A^3*c)+ 2* (b'*A^4*c);
    coneq(end+1) = 6* (3* (b'*A*c)+ ((b'.*c')*A*c)+ 3* ( (b'*A^2*c.^2)+ 4* (b'*diag(c)*A^3*c))) - 2 - 3* (b'*A*c.^2)- 24* (b'*A^2*c)- 36* ((b.*c)'*A^2*c)- 36* (b'*A^3*c.^2);
    coneq(end+1) = 1 - 12* (b'*A*c)+ 24* ((b'.*c')*A*c)- 36* ((b.*c.^2)'*A*c)+ 24* (b'*A*c.^2)- 12* (b'*A*c.^3)+  24* (b'*A^2*c)- 72* ((b.*c)'*A^2*c)+ 72* (b'*diag(c).^2*A^2*c)- 36* (b'*A^2*c.^2)+ 24* (b'*A^2*c.^3);
    coneq(end+1) = 12 *(2* ((b'.*c')*A*c)+ 6* ((b.*c.^2)'*A*c)+ 6* (b'*(A*c).^2)+ 5* (b'*A*c.^2)+ 2* (b'*A^2*c)+ 12* b'*A*diag(c).^2*A*c) - 5 - 24* (b'*A*c)- 144* (b'*diag(c)*(A*c).^2)- 72* (b'*A*c.^3)- 144* (b'*A*diag(c)*A*c);
    coneq(end+1) = 1 + 60* ((b'.*c')*A*c)+ 120* (b'*diag(c).^3*A*c)+ 60* (b'*A*c.^3)- 30* (6* ((b.*c.^2)'*A*c)+ (b'*A*c.^2)+ (b'*A*c.^4));
    coneq(end+1) = 1 + 12* ((b'.*c')*A*c)+ 36* (b'*diag(c).^2*A*c.^2)+ 12* (b'*A*c.^3)-6* (6* ((b.*c.^2)'*A*c)+ (b'*A*c.^2)+ 4* (b'*diag(c)*A*c.^3));
    coneq(end+1) = -1 + 6 * (b'*A*c)+ 6* ((b'.*c')*A*c)- 3* (b'*A*c.^2)+ 12* (b'*A^2*c)- 36* ((b.*c)'*A^2*c)- 18* (b'*A^2*c.^2)+ 36* (b'*diag(c)*A^2*c.^2);
    coneq(end+1) = -(5/144) - (b'*A*c)/6 + (2* ((b'.*c')*A*c))/3 - (b'*(A*c).^2)/2 + (b'*A*c.^2)/6 + (b'*A^2*c)/6 - (b'*A*diag(c)*A*c)+ (b'*A*(A*c).^2);
    coneq(end+1) = c'.^5*b-1/6;
    coneq(end+1) = 8 *( (b'*A*c)+ 3* ( (b'*(A*c).^2)+ (b'*A^2*c.^2)+ 2* (b'*A*diag(c)*A^2*c))) -  1 - 12* (b'*A*c.^2)- 8 * (b'*A^2*c)- 48* (b'*((A^2*c).*(A*c)))-48* (b'*A^2*diag(c)*A*c);
    coneq(end+1) = 2* (6* ((b'.*c')*A*c)+ 12* (b'*(A*c).^2)+ 9* (b'*A*c.^2)+ 8* (b'*A^2*c)+ 24* (b'*diag(c)*A*diag(c)*A*c)+ 12* (b'*A*diag(c)*A*c.^2)) -1 - 4* (b'*A*c)- 24* ((b.*c)'*A*c.^2)- 24* (b'*((A*c.^2).*(A*c)))- 24* ((b.*c)'*A^2*c)- 48* (b'*A*diag(c)*A*c)- 12* (b'*A^2*c.^2);
end