function multi_rk = multi_stab_fun_2_rk_opt(input_file_name,class,objective,varargin)
%function multi_rk = multi_stab_fun_2_rk_opt(input_file_name,class,objective,varargin)
%
%
% This function calls the rk_opt.m function inside a loop for optimizing
% multiple RK methods given their stability function coefficients.
%
%
% input_file_name: name of the file where the coefficients of the stability 
%                  polynomial function can be found. This file must have a
%                  specific header. For instance:
%   
%                  #stability poly.
%                  xxx
%
%                  #stage	 order	 free params. h	 h/s  iter poly. coeffs.
%
% varargin: can contains the input parameters that will override the 
%           default parameters values in rk_opt. See rk_opt function for 
%           more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load and read the file containing the stability polynomial coefficients
[read_fid,number_stab_poly]= read_header(input_file_name);

% Loop over the stability polynomial
for i_stab_poly = 1:number_stab_poly

    % Read number of stage and order of accuracy of the RK scheme
    [s,p,fp,d]= read_data(read_fid);
        
    % Set tall tree numbers (indices) and tall tree values that will be 
    % used to enforce the stability coefficients. Indeed, when "s" is 
    % different from "p" the linear stability region is determined by the 
    % value of the additional tall trees.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Check if the minimum number of stages is consistent with the order of
    % the explicit Runge-Kutta method.
    if (strcmp(class,'erk') | strcmp(class(1:2),'2S') | strcmp(class(1:2),'3S'))
        if ((p <= 4) & (s-p)<0)
            disp('Exit from the code.')
            disp('The number of stages must be at least equal to the order of the Runge-Kutta scheme!')
            break;
        elseif ((p>4) & (p<=6) & (s-p)<1)
            disp('Exit from the code.')
            disp('The number of stages must be at least equal to the order of the Runge-Kutta scheme plus one!')
            break;
        end
    end
    
    % Polynomial coefficients indices
    poly_coeff_ind = p+1:s;

    % Polynomial coefficients values
    poly_coeff_val = d(7+p+1:length(d));

    % Call to rk_opt
    rk = rk_opt(s,p,class,objective,'poly_coeff_ind',poly_coeff_ind,...
                'poly_coeff_val',poly_coeff_val,varargin{:})
    
end

multi_rk = 1;

end

% =========================================================================


% =========================================================================

function [read_fid,number_stab_poly]= read_header(input_file_name)
%function [read_fid,number_stab_poly]= read_data(input_file_name)
%
% Read name of the variables saved in the file

read_fid = fopen(input_file_name,'r');

% Header
tline = fgets(read_fid); 

% Read number of stability polynomial written in the file 
tline = fgets(read_fid);

% Convert string to floating numbers
d = sscanf(tline,'%f');
number_stab_poly = sscanf(tline,'%f');
number_stab_poly = floor(number_stab_poly);

% Read white line
tline = fgets(read_fid);
 
% Header
tline = fgets(read_fid);

% Read white line
tline = fgets(read_fid);

end
% =========================================================================


% =========================================================================

function [s,p,fp,d]= read_data(read_fid)
%function [s,p,fp]= read_data(read_fid)
%
% Read number of stage and order of accuracy of the method

% Read information
tline = fgets(read_fid);
    
% Convert string to floating numbers
d = sscanf(tline,'%f');
    
% Number of stages:
s = floor(d(1));
    
% Order of accuracy:
p = floor(d(2));
    
% Free parameters in the stability function, i.e. difference between the
% actual number of stages and a minimum number of stage which depends on 
% the order of accuracy of the scheme for nonlinear problems.
fp = floor(d(3));

end
% =========================================================================

