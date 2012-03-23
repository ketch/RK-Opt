function multi_rk = multi_rk_opt(inputFileName,class,objective,varargin)
%function multi_rk = multi_rk_opt(inputFileName,class,objective,varargin)
%
%
% This function calls the rk_opt.m function inside a loop for optimizing
% multiple RK methods give their stability function coefficients.
%
%
%
% inputFileName: name of the file where the coefficients of the stability 
%                polynomial function can be found. This file must have a
%                specific header. For instance:
%   
%                #stability poly.
%                18
%
%                #stage	 order	 free params. h	 h/s  iter poly. coeffs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% =========================================================================
% Load and read the file containing the stability polynomial coefficients
% =========================================================================
readFid = fopen(inputFileName,'r');
inputFileName;

% Header
tline = fgets(readFid); 

% Read number of stability polynomial written in the file 
tline = fgets(readFid);

% Convert string to floating numbers
d = sscanf(tline,'%f');
nbrStabPoly = sscanf(tline,'%f');
nbrStabPoly = floor(nbrStabPoly);

% Read White line
tline = fgets(readFid);
 
% Header
tline = fgets(readFid);

% Read White line
tline = fgets(readFid);


% Loop over the stability polynomial
for i_stabPoly = 1:nbrStabPoly
    
    % Read information
    tline = fgets(readFid);
    
    % Convert string to floating numbers
    d = sscanf(tline,'%f');
    
    % Number of stages:
    s = floor(d(1));
    
    % Order of accuracy:
    p = floor(d(2));
    
    % Free parameters in the stability function, i.e. difference between 
    % number of stages (which must be larger than a minimum value) and
    % order of accuracy of the scheme.
    fp = floor(d(3));
    
    % Set tall tree numbers (indices) and tall tree values that will be 
    % used to enforce the stability coefficients. Indeed, when "s" is 
    % different from "p" the linear stability region is determined by the 
    % value of the additional tall trees.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Check if the minimum number of stages is consistent with the order of
    % the Runge-Kutta method.
    if strcmp(class,'erk') 
        if ((p <= 4) & (s-p)<0)
            disp('Exit from the code.')
            disp('The number of stages must be at least equal to the order of the Runge-Kutta scheme!')
            break;
        elseif ((p>4) & (p<=6) & (s-p)<1)
            disp('Exit from the code.')
            disp('The number of stages must be at least equal to the order of the Runge-Kutta scheme plus one!')
            break;
        end
        
        poly_coeff_ind = p+1:s
        poly_coeff_val = d(7+p+1:length(d))
    end
            
    rk = rk_opt(s,p,class,objective,'poly_coeff_ind',poly_coeff_ind,'poly_coeff_val',poly_coeff_val,varargin{:})
    
end

multi_rk = 1;








