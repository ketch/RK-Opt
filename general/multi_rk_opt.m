function multi_rk = multi_rk_opt(inputFileName,class,objective,startvec,solveorderconditions,np,max_tries,writeToFile)
%function multi_rk = multi_rk_opt(inputFileName,class,objective,startvec,solveorderconditions,np,max_tries,writeToFile)
%
%
% This function calls the rk_opt.m function inside a loop for optimizing
% multiple RK methods.
%
%
%
% inputFileName: name of the file where the coefficients of the stability 
%                polynomial function can be found
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
    
    % Free parameters in the stability function optimization
    fp = floor(d(3));
    
    % Set tall tree numbers (indices) and tall tree values
    if fp == 0
        poly_coeff_ind = []
        poly_coeff_val = []
    else
        poly_coeff_ind = s-fp+1:s
        poly_coeff_val = d(6+s-fp+2:length(d))
    end

    
    rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,np,max_tries,writeToFile)
    
end

multi_rk = 1;









