function [h,poly_coeff] = opt_poly_bisect(lam,s,p,basis,varargin)
% function [h,poly_coeff] = opt_poly_bisect(lam,s,p,basis,varargin)
%
% Finds an optimally stable polynomial of degree s and order p for the spectrum
% lam in the interval (h_min,h_max) to precision eps.

[lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,...
        h_true] = opt_poly_params(s,lam,varargin);

diagnostics_init;

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
        lam = lam_func(h)/h;
    end
    
    [status, a, v, diag_solve] = least_deviation(h,lam,s,p,basis,'sdpt3',row_scale,diag_on);
    if strcmp(status,'Failed') || v==Inf || isnan(v)
        disp('sdpt3 failed!');
        [status, a, v, diag_solve] = least_deviation(h,lam,s,p,basis,'sedumi',row_scale,diag_on);
        if strcmp(status,'Failed') || v==Inf || isnan(v)
            disp('sedumi failed!'); 
        end
    end

    diagnostics_update;

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

diagnostics_finalize;

h=h_min; % Return the largest known feasible value



%====================================================
function [status,poly_coeffs,v,diag] = ...
            least_deviation(h,lam,s,p,basis,solver,row_scaling,diag_on,precision)
% function [status,poly_coeffs,v,diag] = ...
%            least_deviation(h,lam,s,p,basis,solver,row_scaling,diag_on,precision)
%
% Solve the least deviation problem \min |R(h\lambda)|
%   for specified set of \lambda, with order constraints: R(z)\approx \exp(z)
%
% Note that poly_coeffs contains the polynomial coefficients
% in the modified/scaled basis!
% Perhaps we should return monomial basis coefficients instead.

if nargin<9; precision='best'; end
if nargin<8, diag_on=0; end
if nargin<7; row_scaling=ones(1,s+1); end
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

% Here we convert to the monomial basis, for convenience of plotting.
% But for high-degree polynomials, it's numerically better to work in the adapted basis.
if ~strcmp(basis,'monomial')
    poly_coeffs = poly_coeffs'*b;
end
status=cvx_status;
v = cvx_optval;
diagnostics_least_deviation;



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



%==============================================================
function [lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,...
        h_true] = opt_poly_params(s,lam,optional_params);
% function [lam_func,tol_bisect,tol_feasible,h_min,h_max,max_steps,...
%        h_true] = opt_poly_params(s,lam,optional_params);
%
% Set default optional and parameter values
i_p = inputParser;
i_p.FunctionName = 'opt_poly_params';

default_h_max        = 2.01*s^2*max(abs(lam)); 

i_p.addParamValue('lam_func',0);
i_p.addParamValue('tol_bisect',1.e-3,@isnumeric);
i_p.addParamValue('tol_feasible',1.e-9,@isnumeric);
i_p.addParamValue('h_min',0,@isnumeric);
i_p.addParamValue('h_max',default_h_max,@isnumeric);
i_p.addParamValue('max_steps',1000,@isnumeric);
i_p.addParamValue('h_true',[]);

i_p.parse(optional_params{:});

lam_func          = i_p.Results.lam_func;
tol_bisect        = i_p.Results.tol_bisect;
tol_feasible      = i_p.Results.tol_feasible;
h_min             = i_p.Results.h_min;
h_max             = i_p.Results.h_max;
max_steps         = i_p.Results.max_steps;
h_true            = i_p.Results.h_true;
