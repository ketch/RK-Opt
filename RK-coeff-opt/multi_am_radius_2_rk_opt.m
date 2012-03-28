function multi_rk = multi_am_radius_2_rk_opt(input_file_name,class,varargin)
%function multi_rk = multi_am_radius_2_rk_opt(input_file_name,class,varargin)
%
% This function calls the rk_opt.m function inside a loop for optimizing
% multiple explicit monotonicity preserving RK methods given their stability 
% function coefficients. 
%
% Currently the stability function coefficients in rk_opt.m is imposed using the
% monomial basis, i.e. \Psi(z) = sum_{j=0}^m a_j z^j.
%
%
% input_file_name: name of the file where the coefficients of the stability 
%                  polynomial function can be found. This file must have a
%                  specific structure. For instance:
%   
%                  stability polynomials
%                  1
%
%                  ==========================================================
%
%                  order
%                  2.0000000000000000E+00	
%
%                  step
%                  1.0000000000000000E+00	
%
%                  stage
%                  2.0000000000000000E+00	
%
%                  R
%                  1.0000000099992548E+00	
%
%                  gamma
%                  2.0599050593350069E-15	0.0000000000000000E+00	9.9999999999999789E-01	
%
%                  coeffs_monomial
%                  1.0000000000000000E+00	1.0000000000000000E+00	2.5000000000000056E-01	
%
%                  ==========================================================
%
%
% varargin: can contains the input parameters that will override the default 
%           parameters values in rk_opt. See rk_opt function for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up parameters
[problem_class]= setup_params(varargin);

% Load and read the file containing the stability polynomial coefficients
[read_fid,number_stab_poly]= read_header(input_file_name);


% Loop over the stability polynomial
for i_stab_poly = 1:number_stab_poly
    
    % Read data from file
    [s,k,p,coeffs_monomial]= read_data(read_fid,class);
   

    % Call to rk_opt
    if strcmp(problem_class,'linear')
        
        % The order of accuracy is set to one to re-use completely the 
        % functions implemented in the RK-coeff-opt to impose just the 
        % linear order conditions.
        h_p = 1
        
        % Polynomial coefficients indicies
        poly_coeff_ind = 2:s;
        
        % Polynomial coefficients value
        poly_coeff_val = coeffs_monomial(3:length(coeffs_monomial));

    else
        h_p = p;

        % Polynomial coefficients indicies
        poly_coeff_ind = p+1:s;
        
        % Polynomial coefficients value
        poly_coeff_val = coeffs_monomial(p+2:length(coeffs_monomial));

    end

    objective = 'acc'; % Optimize the leading truncation error constant.
    rk = rk_opt(s,h_p,class,objective,'poly_coeff_ind',poly_coeff_ind,'poly_coeff_val',poly_coeff_val,varargin{:})
    
end

multi_rk = 1;

end
% =========================================================================


% =========================================================================

function [problem_class]= setup_params(optional_params)
%function [problem_class]= setup_params(optional_params)
%
% Set default optional and param values

i_p = inputParser;
i_p.FunctionName = 'setup_params';

% Expected values
expected_problem_class = ['linear','nonlinear'];

% Default values
default_problem_class = 'linear';


% Populate input parser object
% ----------------------------
% Parameter values
i_p.addParamValue('problem_class',default_problem_class,@(x) iscahr(x) && any(validatestring(x,expected_problem_class))); 

i_p.parse(optional_params{:});

% Dereference parameters
problem_class = i_p.Results.problem_class;

end
% =========================================================================


% =========================================================================

function [read_fid,number_stab_poly]= read_header(input_file_name)
%function [read_fid,number_stab_poly]= read_data(input_file_name)
%
% Read the number of stability polynomials saved in the file

read_fid = fopen(input_file_name,'r');

% Header
tline = fgets(read_fid); 

% Read number of stability polynomial written in the file 
tline = fgets(read_fid);

% Convert string to floating numbers
d = sscanf(tline,'%f');
number_stab_poly = sscanf(tline,'%f');
number_stab_poly = floor(number_stab_poly);


end
% =========================================================================


% =========================================================================

function [s,k,p,coeffs_monomial]= read_data(read_fid,class)
%function [s,k,p,coeffs_monomial]= read_data(read_fid,class)
%
% Read number of stage and order of the method

    % Read separation line
    tline = fgets(read_fid);    
    
    % Read white line
    tline = fgets(read_fid);

    % Read string 'order'
    tline = fgets(read_fid);

    % Read value of the order
    tline = fgets(read_fid);
    
    % Convert string to floating numbers
    p = floor(sscanf(tline,'%f'));
    
    % Read white line
    tline = fgets(read_fid);

    % Read string 'step'
    tline = fgets(read_fid);

    % Read # of the step
    tline = fgets(read_fid);

    % Convert string to floating numbers
    k = floor(sscanf(tline,'%f'));

    % Read white line
    tline = fgets(read_fid);

    % Read string 'stage'
    tline = fgets(read_fid);

    % Read # of the stages
    tline = fgets(read_fid);

    % Convert string to floating numbers
    s = floor(sscanf(tline,'%f'));

    % Check if the minimum number of stages is consistent with the order of
    % the explicit Runge-Kutta method.
    if (strcmp(class,'erk') | strcmp(class(1:2),'2S') | strcmp(class(1:2),'3S'))
        if ((p <= 4) & (s-p)<0)
            disp('Exit from the code.')
            disp('The number of stages must be at least equal to the order of the Runge-Kutta scheme!')
            return;
        elseif ((p>4) & (p<=6) & (s-p)<1)
            disp('Exit from the code.')
            disp('The number of stages must be at least equal to the order of the Runge-Kutta scheme plus one!')
            return;
        end
    end


    % Read white line
    tline = fgets(read_fid);

    % Read string 'R'
    tline = fgets(read_fid);

    % Read value of R
    tline = fgets(read_fid);

    % Convert string to floating numbers
    R = floor(sscanf(tline,'%f'));

    % Read white line
    tline = fgets(read_fid);

    % Read string 'gamma'
    tline = fgets(read_fid);

    % Read components of the array gamma
    tline = fgets(read_fid);

    % Convert string to floating numbers
    d = sscanf(tline,'%f');
    gamma = d(1:length(d));

    % Read white line
    tline = fgets(read_fid);

    % Read string 'coeffs_monomial'
    tline = fgets(read_fid);

    % Read components of the array coeffs_monomial
    tline = fgets(read_fid);

    % Convert string to floating numbers
    d = sscanf(tline,'%f');
    coeffs_monomial = d(1:length(d));

end
% =========================================================================
