function [h,poly_coeff,diag] = opt_poly_bisect(lam,s,p,basis,lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,initial_scaling,h_true)
%function [h,poly_coeff,diag] = opt_poly_bisect(lam,s,p,basis,lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,initial_scaling,h_true)
%
% Finds an optimally stable polynomial of degree s and order p for the spectrum
% lam in the interval (h_min,h_max) to precision eps.

if nargin<5 lam_func=0; end
if nargin<6 basis='monomial'; end
if nargin<7 tol_bisect=1.e-2; end
if nargin<8 tol_feasible=1.e-10; end
if nargin<9 h_min=0; end
if nargin<10 h_max=100; end
if nargin<11 max_steps=1000; end
if nargin<12 
    if strcmp(basis,'monomial')
        initial_scaling = factorial(0:s);
    else
        initial_scaling = ones(s+1,1);
    end
end
if nargin<13 h_true=[]; end

if nargout >3           %Diagnostics requested
    if nargin < 13
        error('Enabling diagnostics assumes h_true is given!');
    end
    diag_on=true;
    diag.lam = lam; diag.s = s; diag.p = p; diag.basis = basis; diag.tol_bisect = tol_bisect;
    diag.tol_feasible = tol_feasible; diag.h_min = h_min; diag.h_max = h_max; diag.max_steps = max_steps; diag.h_true = h_true;
    diag.solved_steps = []; diag.failed_steps = []; diag.right_steps = []; diag.wrong_steps = [];
else
    diag_on=false;
end

scale=initial_scaling;
for i_step=1:max_steps
    %Stop bisecting when tolerance is achieved
    if ((h_max-h_min)/h_min < tol_bisect) || (h_max < tol_bisect)
        break;
    end

    h = (h_max + h_min)/2.;

    if isfloat(lam_func) == 0
        lam = lam_func(h)/h;
    end
    
    if strcmp(basis,'monomial')
        [status, x_opt, v, diag_solve] = monomial_solve(h,lam,s,p,'sdpt3',scale);
    elseif strcmp(basis,'modified')
        [status, x_opt, v, diag_solve] = modified_basis_solve(h,lam,s,p,basis,'sdpt3');
        if strcmp(status,'Failed') || v==Inf
            disp('sdpt3 failed!');
            [status, x_opt, v, diag_solve] = modified_basis_solve(h,lam,s,p,basis,'sedumi');
            if strcmp(status,'Failed') || v==NaN
                disp('sedumi failed!'); 
            end
        end
    end

    if diag_on
        diag.step(i_step).solve = diag_solve; diag.step(i_step).h = h; diag.step(i_step).h_max = h_max;
        diag.step(i_step).h_min = h_min; diag.step(i_step).x_opt = x_opt; diag.step(i_step).v = v;
        diag.step(i_step).status = status;
    end

    % Check if CVX has found a feasible solution and bisect accordingly
    if strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved')
        if v<=tol_feasible
            h_min = h;
            if strcmp(basis,'monomial')
                poly_coeff = [1./factorial(0:p) x_opt'./scale(p+2:end)]';
                scale(p+2:end)=1./poly_coeff';
                scale(p+2:end)=scale(p+2:end)./x_opt;
            else
                poly_coeff = x_opt;
            end
        else %Infeasible
            h_max = h;
        end
        
        if diag_on
            diag.solved_steps(end+1) = i_step;
        end

    else % Solver Failed; TODO: set some status flag
        h_max = h;
        if diag_on
            diag.failed_steps(end+1) = i_step;
        end
    end
    
    if diag_on
        if h_true > h_min && h_true < h_max
            diag.right_steps(end+1) = i_step;
        else
            diag.wrong_steps(end+1) = i_step;
        end
    end
    fprintf(' [%5.0d] h_min: %e h_max: %e h_true: %e\n',i_step,h_min,h_max,h_true);         
    h/s^2
end

h=h_min; % Return the largest known feasible value

if diag_on
    diag.final_solve.solve = diag_solve; diag.final_solve.h = h; diag.final_solve.h_max = h_max;
    diag.final_solve.h_min = h_min; diag.final_solve.x_opt = x_opt; diag.final_solve.v = v;
    diag.final_solve.status = status;
    diag.n_iter = i_step;
end

%====================================================
function [status,poly_coeffs,v,diag] = ...
            least_deviation_modified_basis(h,lam,s,p,basis,solver,scaling,precision)
%function [status,poly_coeffs,v,diag] = ...
%            least_deviation_modified_basis(h,lam,s,p,basis,solver,scaling,precision)
%
% Solve the least deviation problem \min |R(h\lambda)|
%
% Note that poly_coeffs contains the polynomial coefficients
% in the modified basis!

if nargin<8 precision='best'; end
if nargin<7 scaling=ones(1,s+1); end
if nargin<6 solver = 'sdpt3'; end
if nargin<5 basis = 'chebyshev'; end

min_real_part = min(real(h*lam));
max_imag_part = max(abs(imag(h*lam)));
if strcmp(basis,'chebyshev')
    [b,c] = scaled_chebyshev_basis(s,min_real_part,0,h*lam);
elseif strcmp(basis,'binomial')
    [b,c] = scaled_binomial_basis(s,-min_real_part/2,h*lam);
elseif strcmp(basis,'rotated chebyshev')
    [b,c] = scaled_chebyshev_basis(s,-1i*max_imag_part,1i*max_imag_part,h*lam);
    %Now rotate (make all basis coefficients real):
    for i=1:(s+1)/2
      b(2*i,:)=1i*b(2*i,:);
      c(:,2*i) = 1i*c(:,2*i);
    end
end

order_conditions = scaling(1:p+1)'./factorial(0:p)';

for i=1:s+1
    b(:,i) = b(:,i)./scaling';
    c(i,:) = c(i,:)./scaling;
end
        
cvx_begin
  cvx_quiet(true)
  cvx_precision(precision)
  cvx_solver(solver)
  
  variable poly_coeffs(s+1) 
  
  R=abs(c*poly_coeffs)-1;
  
  b(:,1:p+1)'*poly_coeffs==order_conditions;
  minimize max(R)
cvx_end

status=cvx_status;
v = cvx_optval;

if nargout>3
    diag.b = b;
    diag.c = c;
    diag.order_conditions = order_conditions;
    diag.cond_b = cond(b);
    diag.cond_c = cond(c);
    diag.h = h;
    diag.precision = precision;
    diag.poly_coeffs = poly_coeffs;
    diag.cvx.cputime = cvx_cputime;
    diag.cvx.optdpt = cvx_optdpt;
    diag.cvx.optpnt = cvx_optpnt;
    diag.cvx.slvitr = cvx_slvitr;
    diag.cvx.slvtol = cvx_slvtol;
    diag.cvx.status = cvx_status;
end
