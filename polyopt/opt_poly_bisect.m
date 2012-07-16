function [h,poly_coeff] = opt_poly_bisect(lam,s,p,basis,lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,initial_scaling,h_true)
%function [h,poly_coeff] = opt_poly_bisect(lam,s,p,basis,lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,initial_scaling,h_true)
%
% Finds an optimally stable polynomial of degree s and order p for the spectrum
% lam in the interval (h_min,h_max) to precision eps.

if nargin<4; basis='monomial'; end
if nargin<5; lam_func=0; end
if nargin<6; tol_bisect=1.e-3; end
if nargin<7; tol_feasible=1.e-9; end
if nargin<8; h_min=0; end
if nargin<9  h_max=2.01*s^2*max(abs(lam)); end
if nargin<10; max_steps=1000; end
if nargin<11 
    if strcmp(basis,'monomial')
        initial_scaling = factorial(0:s);
    else
        initial_scaling = ones(1,s+1);
    end
end
if nargin<12; h_true=[]; end

diagnostics_init;

scale=initial_scaling;
for i_step=1:max_steps
    %Stop bisecting when relative tolerance is achieved
    if ((h_max-h_min)/h_min < tol_bisect) || (h_max < tol_bisect)
        break;
    end

    h = (h_max + h_min)/2.;

    if isfloat(lam_func) == 0
        lam = lam_func(h)/h;
    end
    
    [status, a, v, diag_solve] = least_deviation(h,lam,s,p,basis,'sdpt3',diag_on,scale);
    if strcmp(status,'Failed') || v==Inf
        disp('sdpt3 failed!');
        [status, a, v, diag_solve] = least_deviation(h,lam,s,p,basis,'sedumi',diag_on,scale);
        if strcmp(status,'Failed') || isnan(v)
            disp('sedumi failed!'); 
        end
    end

    diagnostics_update;

    % Check if CVX has found a feasible solution and bisect accordingly
    if strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved')
        if v<=tol_feasible
            h_min = h;
            if strcmp(basis,'monomial')
                poly_coeff = [1./factorial(0:p) a'./scale(p+2:end)];
                scale(p+2:end)=scale(p+2:end)./a';
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

diagnostics_finalize;

h=h_min; % Return the largest known feasible value



%====================================================
function [status,poly_coeffs,v,diag] = ...
            least_deviation(h,lam,s,p,basis,solver,diag_on,row_scaling,precision)
% function [status,poly_coeffs,v,diag] = ...
%            least_deviation(h,lam,s,p,basis,solver,diag_on,row_scaling,precision,diag_on)
%
% Solve the least deviation problem \min |R(h\lambda)|
%
% Note that poly_coeffs contains the polynomial coefficients
% in the modified/scaled basis!
% Perhaps we should return monomial basis coefficients instead.

if nargin<9; precision='best'; end
if nargin<8; row_scaling=ones(1,s+1); end
if nargin<7, diag_on=0; end
if nargin<6; solver = 'sdpt3'; end
if nargin<5; basis = 'chebyshev'; end

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
  
    if strcmp(basis,'monomial')
        variable poly_coeffs(s-p);
        fixedvec = c(:,1:p+1)*fixed_coefficients;
        R=abs(fixedvec+c(:,p+2:end)*poly_coeffs)-1.;
    else
        variable poly_coeffs(s+1) 
        b(:,1:p+1)'*poly_coeffs==fixed_coefficients;
        R=abs(c*poly_coeffs)-1;
    end
    minimize max(R)
cvx_end

status=cvx_status;
v = cvx_optval;
diagnostics_least_deviation;
