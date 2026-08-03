function [con,coneq]=nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k,emb_poly_coeff_ind,emb_poly_coeff_val,constrain_emb_stability,c_lower_bound,c_upper_bound,c_monotone,interp_nodes,interp_values,interp_eval)
% function [con,coneq]=nonlinear_constraints(x,class,s,p,objective,poly_coeff_ind,poly_coeff_val,k,emb_poly_coeff_ind,emb_poly_coeff_val,constrain_emb_stability,c_lower_bound,c_upper_bound,c_monotone,interp_nodes,interp_values,interp_eval)
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
%     * *c_lower_bound*: lower bound of the Butcher coefficients c
%     * *c_upper_bound*: upper bound of the Butcher coefficients c
%     * *c_monotone*: set this to true if the Butcher coefficients c should be monotonically increasing
%     * *interp_nodes*: nodes used for stability-function matching.
%     * *interp_values*: values of the stability function at interp_nodes
%     * *interp_eval*: evaluation used for the interpolation constraints, either 'butcher' or 'lowstorage'.
%
% The outputs are:
%     * *con*: inequality constraints, i.e. absolute monotonicity conditions if objective = 'ssp' or nothing if objective = 'acc'
%     * *coneq*: order conditions plus stability function coefficients constraints (tall-tree elementary weights) and interpolation constraints
% 
% Two forms of the order conditions are implemented: one based on **Butcher's
% approach**, and one based on **Albrecht's approach**. One or the other may lead 
% to a more tractable optimization problem in some cases, but this has not been 
% explored carefully. The Albrecht order conditions are implemented up to order 9, assuming
% a certain stage order, while the Butcher order conditions are implemented up to order 9 but
% do not assume anything about the stage order. Albrecht's approach is used
% by default.


oc_form = 'albrecht';

% If interp_eval is not provided, use the classic Butcher-form evaluation by default.
if nargin < 17 || isempty(interp_eval)
    interp_eval = 'butcher';
end

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
% Interpolation constraints for the stability function
if ~isempty(interp_nodes)
    z_nodes  = interp_nodes(:);
    R_target = interp_values(:);

    if numel(z_nodes) ~= numel(R_target)
        error('interp_nodes and interp_values must have same length.');
    end

    R_rk = zeros(numel(z_nodes),1);

    switch lower(interp_eval)
    case 'butcher'
        % Evaluate R(z) from the Butcher coefficients.
        e = ones(s,1);
        for j = 1:numel(z_nodes)
            z = z_nodes(j);
            Y = (eye(s) - z*A)\e;
            R_rk(j) = 1 + z*(b.'*Y);
        end

    case 'lowstorage'
        if ~strcmp(class,'3Sstar')
            error(['interp_eval = ''lowstorage'' is currently supported ', ...
                   'only for class = ''3Sstar''.']);
        end

        % Evaluate R(z) by applying the 3S* low-storage recursion to the
        % scalar test equation u' = z*u with initial value u^n = 1.
        [~,~,~,~,~,~,~,beta,gamma1,gamma2,gamma3,delta] = unpack_lsrk(x,class);
        for j = 1:numel(z_nodes)
            z = z_nodes(j);
            S1 = 1; S2 = 0; S3 = 1;
            for i = 2:(s+1)
                S2 = S2 + delta(i-1)*S1;
                S1 = gamma1(i)*S1 + gamma2(i)*S2 + gamma3(i)*S3 + beta(i,i-1)*z*S1;
            end
            R_rk(j) = S1;
        end

    otherwise
        error('Unknown interp_eval value: %s.', interp_eval);
    end

    % Match the real parts of the stability-function values at the interpolation nodes.
    coneq_interp = (real(R_rk) - real(R_target)).';
    coneq        = [coneq  coneq_interp];
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
%=====================================================
% Constraints on the Butcher parameter c

if ~isempty(c_lower_bound)
    % The solver attempts to satisfy con <= 0 for all entries of con.
    con = [con; c_lower_bound - c];
end

if ~isempty(c_upper_bound)
    % The solver attempts to satisfy con <= 0 for all entries of con.
    con = [con; c - c_upper_bound];
end

if c_monotone
    % The solver attempts to satisfy con <= 0 for all entries of con.
    con = [con; c(1:end-1) - c(2:end)];
end
