function [h,poly_coeff,diag_bisect] = opt_poly_bisect(lam,s,p,basis,varargin)
% function [h,poly_coeff,diag_bisect] = opt_poly_bisect(lam,s,p,basis,varargin)
%
% Finds an optimally stable polynomial of degree s and order p for the spectrum
% lam in the interval (h_min,h_max) to precision eps.
%
% Optional arguments:
%
%
%       solvers: 
%                 A cell array of cvx solver names that should be used to
%                 solve the convex problem in the inner loop. Defaults to
%                 {'sdpt3', 'sedumi'}. You can also add 'mosek' and
%                 'gurobi' if you have obtained (free academic) licences
%                 for these and installed them in cvx.
%       lam_func: 
%                 A function used to generate the appropriate spectrum
%                 at each bisection step, instead of using a fixed (scaled) spectrum.
%                 Used for instance to find the longest rectangle of a fixed height
%                 (see Figure 10 of :cite:`2012_optimal_stability_polynomials`).
%       eta:
%                 A damping parameter, which forces the polynomial to
%                 satisfy |R(h*lambda)| < 1-eta instead of '... < 1'
%                 This ensures, that there is a small stable environment
%                 around each eigenvalue
%
%
% Examples:
%
%       - To find negative real axis inclusion::
%
%               lam = spectrum('realaxis',500);       
%               s = 10; p = 2;
%               [h,poly_coeff] = opt_poly_bisect(lam,s,p,'chebyshev')    
%
%       - To reproduce figure 10 of :cite:`2012_optimal_stability_polynomials` ::
%
%               lam_func = @(kappa) spectrum('rectangle',100,kappa,10)
%               [h,poly_coeff] = opt_poly_bisect(lam,20,1,'chebyshev','lam_func',lam_func)
%               plotstabreg_func(poly_coeff,[1])

[lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,...
        h_true,do_plot,solvers,eta] = opt_poly_params(s,lam,varargin);

%Diagnostics requested
if nargout >= 3           
    if isempty(h_true)
        error('Enabling diagnostics assumes h_true is given!');
    end
    diag_on=true;
    diag_bisect.lam = lam; diag_bisect.s = s; diag_bisect.p = p; diag_bisect.basis = basis; diag_bisect.tol_bisect = tol_bisect;
    diag_bisect.tol_feasible = tol_feasible; diag_bisect.h_min = h_min; diag_bisect.h_max = h_max; diag_bisect.max_steps = max_steps; diag_bisect.h_true = h_true;
    diag_bisect.solved_steps = []; diag_bisect.failed_steps = []; diag_bisect.right_steps = []; diag_bisect.wrong_steps = [];
else
    diag_on=false;
end


if strcmp(basis,'monomial')
    row_scale = factorial(0:s);
else
    row_scale = ones(1,s+1);
end

for i_step=1:max_steps
    %Stop bisecting when relative tolerance is achieved
    if ((h_max-h_min)/h_min < tol_bisect) || (h_max < tol_bisect)
        break;
    end

    h = (h_max + h_min)/2.;

    if isfloat(lam_func) == 0
        % Divide by 'h' here, since we will multiply by it later
        lam = lam_func(h)/h;
    end

    % lam should be a column vector
    shape = size(lam);
    if shape(1) == 1
        lam = lam';
    end
    % length(lam) should be at least s+1
    assert(length(lam)>=s+1,'Underdetermined: spectrum should contain at least s+1 values.')

    
    for solver_idx = 1:length(solvers)
        solver = solvers{solver_idx};  
        [status, a, v, diag_solve] = least_deviation(h,lam,s,p,eta,basis,solver,row_scale,diag_on);
        
        if strcmp(status,'Failed') || v==Inf || isnan(v)
            fprintf('%s failed!\n', solver);
        else
            break
        end
    end

    if diag_on
        diag.step(i_step).solve = diag_solve; diag.step(i_step).h = h; diag.step(i_step).h_max = h_max;
        diag.step(i_step).h_min = h_min; diag.step(i_step).x_opt = x_opt; diag.step(i_step).v = v;
        diag.step(i_step).status = status;

        if strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved')
            diag.failed_steps(end+1) = i_step;
        end
 
    end

    % Check if CVX has found a feasible solution and bisect accordingly
    if strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved')
        if v<=tol_feasible
            h_min = h;
            if strcmp(basis,'monomial')
                poly_coeff = [1./factorial(0:p) a'./row_scale(p+2:end)];
                row_scale(p+2:end)=row_scale(p+2:end)./a';
            else
                poly_coeff = a;
            end
        else %Infeasible
            h_max = h;
        end
    else % Solver Failed; TODO: set some status flag
        h_max = h;
    end
    
    fprintf(' [%5.0d] h_min: %e h_max: %e h_true: %e\n',i_step,h_min,h_max,h_true);         
end

if diag_on
    diag.final_solve.solve = diag_solve; diag.final_solve.h = h; diag.final_solve.h_max = h_max;
    diag.final_solve.h_min = h_min; diag.final_solve.x_opt = x_opt; diag.final_solve.v = v;
    diag.final_solve.status = status;
    diag.n_iter = i_step;
end

h=h_min; % Return the largest known feasible value

if do_plot
    stability_plot(h*lam,poly_coeff);
end
end


%====================================================
function [status,poly_coeffs,v,diag] = ...
            least_deviation(h,lam,s,p,eta,basis,solver,row_scaling,diag_on,precision)
% function [status,poly_coeffs,v,diag] = ...
%            least_deviation(h,lam,s,p,basis,solver,row_scaling,diag_on,precision)
%
% Solve the least deviation problem \min |R(h\lambda)|
%   for specified set of \lambda, with order constraints: R(z)\approx \exp(z)
%
% Note that poly_coeffs contains the polynomial coefficients
% in the modified/scaled basis!
% Perhaps we should return monomial basis coefficients instead.

if nargin<10; precision='best'; end
if nargin<9, diag_on=0; end
if nargin<8; row_scaling=ones(1,s+1); end
if nargin<7; solver = 'sdpt3'; end
if nargin<6; basis = 'chebyshev'; end

% b(j,:): coefficients of the jth basis polynomials in terms of monomials
% c(i,j): jth basis polynomial evaluated at h*lam
min_real_part = min(real(h*lam));
max_imag_part = max(abs(imag(h*lam)));
if strcmp(basis,'chebyshev')
    assert(min_real_part~=0,'Use chebyshev basis only with spectra that have negative real part.')
    [b,c] = scaled_chebyshev_basis(s,min_real_part,0,h*lam);
elseif strcmp(basis,'binomial')
    [b,c] = scaled_binomial_basis(s,-min_real_part/2,h*lam);
elseif strcmp(basis,'rotated chebyshev')
    assert(max_imag_part~=0,'Use rotated chebyshev basis only with spectra that have nonzero imaginary part.')
    [b,c] = scaled_chebyshev_basis(s,-1i*max_imag_part,1i*max_imag_part,h*lam);
    %Now rotate (make all basis coefficients real):
    for i=1:(s+1)/2
      b(2*i,:)=1i*b(2*i,:);
      c(:,2*i) = 1i*c(:,2*i);
    end
elseif strcmp(basis,'monomial')
    c = zeros(length(lam),s+1);
    for i=1:length(lam)
        c(i,:) = (h*lam(i)).^(0:s)./row_scaling;
    end
end

fixed_coefficients = row_scaling(1:p+1)'./factorial(0:p)';

if ~strcmp(basis,'monomial')
    for i=1:s+1
        c(i,:) = c(i,:)./row_scaling;
        b(:,i) = b(:,i)./row_scaling';
    end
end
   
cvx_begin
    cvx_quiet(true)
    cvx_precision(precision)
    cvx_solver(solver)
  
    max_value = 1.0-eta; % previously 1.0
    if strcmp(basis,'monomial')
        variable poly_coeffs(s-p);
        fixedvec = c(:,1:p+1)*fixed_coefficients;
        R=abs(fixedvec+c(:,p+2:end)*poly_coeffs)-max_value; %%damping
    else
        variable poly_coeffs(s+1) 
        b(:,1:p+1)'*poly_coeffs==fixed_coefficients;
        R=abs(c*poly_coeffs)-max_value;
    end
    minimize max(R)
cvx_end

% Here we convert to the monomial basis, for convenience of plotting.
% But for high-degree polynomials, it's numerically better to work in the adapted basis.
if ~strcmp(basis,'monomial')
    poly_coeffs = poly_coeffs'*b;
end
status=cvx_status;
v = cvx_optval;
if diag_on
    diag.b = b;
    diag.c = c;
    diag.fixed_coefficients = fixed_coefficients;
    diag.cond_b = cond(b);
    diag.cond_c = cond(c);
    diag.h = h;
    diag.precision = precision;
    diag.poly_coeffs = poly_coeffs;
    diag.cvx.cputime = cvx_cputime;
    diag.cvx.slvitr = cvx_slvitr;
    diag.cvx.slvtol = cvx_slvtol;
    diag.cvx.status = cvx_status;
else 
    diag = [];
end
end


%====================================================
function [b,c] = scaled_binomial_basis(N,r,z)
% function [b,c] = scaled_binomial_basis(N,r,z)
% Given an arbitrary radius r, scaled_binomial generates a basis of polynomials (1+z/r)^j. 
% N is the order of polynomial basis desired.  
%
% Returns:
%  1. A matrix b, whose jth row contains the coefficients of the jth basis function
%     with the coefficients in order of ascending degree
%
%  2. A matrix c, whose jth column contains the values of the jth basis function
%     evaluated at the points z.
%
% Some of the loops could be vectorized, but since this routine isn't a bottleneck we've
% opted for clarity.

b = zeros(N+1,N+1);

for j=0:N
    for k=0:j
        b(j+1,k+1) = nchoosek(j,k)/r^k;
    end
end

if nargin < 3
    return
end

c = zeros(length(z),N+1);

for j=0:N
    for i=1:length(z)
        c(i,j+1) = (1+z(i)/r)^j;
    end
end
end


%====================================================
function [b,c] = scaled_chebyshev_basis(N,zmin,zmax,z)
% function [b,c] = scaled_chebyshev_basis(N,zmin,zmax,z)
% Given an arbitrary domain on the real axis [zmin,zmax], scaled_chebyshev generates a basis of Chebyshev Polynomials of
% the first kind scaled and shifted by the affine mapping: m(x)=m1*x+m0 where m1=2/(zmax-zmin) and m0=-(1+zmin) so m([zmin,zmax]) -> [-1,1]
% N is the order of polynomial basis desired.  
%
% Returns:
%  1. A matrix b, whose ith row contains the coefficients of the ith basis function
%     with the coefficients in order of ascending degree
%
%  2. A matrix c, whose ith column contains the values of the ith basis function
%     evaluated at the points z.
%
% Some of the loops could be vectorized, but since this routine isn't a bottleneck we've
% opted for clarity.

b = zeros(N+1,N+1);

m1 = 2/(zmax-zmin);
m0 = -(1+zmin*m1);

% T_0' = 1
b(1,1) = 1;

% T_1' = m1*x + m0
b(2,1) = m0;
b(2,2) = m1;

% T_{n+1}' = 2*(m1*x + m0)*T_{n} - T_{n-1} 
for k=1:N-1 
    b(k+2,:) = 2*(m1*[0 b(k+1,1:end-1)] + m0*b(k+1,:)) - b(k,:);
end

if nargin < 4
    return
end

c = zeros(length(z),N+1);

c(:,1) = 1;
c(:,2) = m1*z+m0;

for k=1:N-1 
    c(:,k+2) = 2*(m1*c(:,k+1).*z + m0*c(:,k+1)) - c(:,k);
end
end


%==============================================================
function [lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,...
        h_true,do_plot,solvers,eta] = opt_poly_params(s,lam,optional_params)
% function [lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,...
%        h_true,do_plot,solvers] = opt_poly_params(s,lam,optional_params);
%
% Set default optional and parameter values
i_p = inputParser;
i_p.FunctionName = 'opt_poly_params';

default_h_max        = 2.01*s^2*max(abs(lam)); 

i_p.addParameter('lam_func',0);
i_p.addParameter('tol_bisect',1.e-3,@isnumeric);
i_p.addParameter('tol_feasible',1.e-9,@isnumeric);
i_p.addParameter('h_min',0,@isnumeric);
i_p.addParameter('h_max',default_h_max,@isnumeric);
i_p.addParameter('max_steps',1000,@isnumeric);
i_p.addParameter('h_true',[]);
i_p.addParameter('do_plot',false,@islogical);
i_p.addParameter('solvers',{'sdpt3','sedumi'});
i_p.addParameter('eta', 0.0);

i_p.parse(optional_params{:});

lam_func          = i_p.Results.lam_func;
tol_bisect        = i_p.Results.tol_bisect;
tol_feasible      = i_p.Results.tol_feasible;
h_min             = i_p.Results.h_min;
h_max             = i_p.Results.h_max;
max_steps         = i_p.Results.max_steps;
h_true            = i_p.Results.h_true;
do_plot           = i_p.Results.do_plot;
solvers           = i_p.Results.solvers;
eta               = i_p.Results.eta;
end


function status = stability_plot(lam,poly_coeff)
plotstabreg_func(poly_coeff,[1])
hold on;
plot(real(lam),imag(lam),'o')
hold off;
status = 1;
end
