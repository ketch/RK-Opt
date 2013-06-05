function [con,coneq]=nonlinear_constraints(x,class,s,p,alpha,beta)
%function [con,coneq]=nonlinear_constraints(x,class,s,p,alpha,beta)
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

[A,b,c]=unpack_rk(x,s,class);

r=-x(end); %Radius of absolute monotonicity

%=====================================================
% Inequality constraints: absolute monotonicity conditions
es=ones(s,1);
K=[A es*0;b' 0];
G=eye(s+1)+r*K;

con1=G\K;
con2=G\[es;1];

con=-[con1(:);con2(:)];
%=====================================================

%=====================================================
% Order conditions
if strcmp(oc_form,'albrecht')
    coneq = oc_albrecht(A,b,c,p);
elseif strcmp(oc_form,'butcher')
    coneq = oc_butcher(A,b,c,p);
end
%=====================================================

% LMM monotonicity constraints
m = 100; % # of multistep method steps to check
k = length(alpha); % # of steps used by multistep method (in one step)
m0 = s+1;         % # of values computed in starting procedure
az = zeros(k,1);
bz = zeros(k,1);
az(1)=alpha(end);
bz(1)=beta(end);
A0 = toeplitz(az,alpha(end:-1:1));
B0 = toeplitz(bz,beta(end:-1:2));
az = zeros(m,1); az(2:k+1) = alpha; AA = toeplitz(az,zeros(1,m));
bz = zeros(m,1); bz(1:k+1) = beta;  bz2 = zeros(1,m); bz2(1)=beta(1); BB = toeplitz(bz,bz2);
J = eye(m); J = J(:,1:k);
R = (eye(m) - AA + r*BB)\J;
J0 = zeros(k,m0);  J0(s,s)=1;
es=ones(s,1);
K0=[A es*0;b' 0];
I = eye(m0);

e0=ones(m0,1);
M1 = R*(A0-r*B0)*J0*((I+r*K0)\e0);
M2 = R*( (A0-r*B0)*J0*r*((I+r*K0)\K0) + r*B0*J0 );
con = [con;-M1(:);-M2(:)];
