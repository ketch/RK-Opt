function multi_R = multi_R_opt(k,p,class,varargin)
%multi_R = multi_R_opt(s,k,p,class,varargin)
%
% This function is an interface to Rskp or Rkp_dw or Rkp_imp and 
% Rkp_imp_dw for running multiple consecutive optimization problems with 
% different values of the input parameters, i.e.: 
%
% k = [k1, k2, ..., kK]^T, K = length(k),  ith-element = # of steps
% p = [p1, p2, ..., pP]^T, P = length(p),  ith-element = order of accuracy
% 
% and 
%
% s = [s1, s2, ..., sS]^T, S = length(s),  ith-element = # of stages
%
% when optimal contractive k-step, s-stage GLM are investigated
%
% The family of method to be considered is specified in the string 'class'.
%
% Note that:
%
% 1- in general S != K != P. Fixed the order of accuracy of the time 
%    integration scheme, one is usually interested in understanding the
%    behavior of the threshold factor R as a function of the number of
%    stages. Therefore, for a fixed element of the array "p", this function
%    loops over the elements of the array "s". Thus, min(s) => max(p). The
%    equality holds for any order of accuracy because the number of 
%    linear order conditions that will be imposed to construct the 
%    GLM coefficients is p. 

% Parse optional input arguments
[s]= setup_params(p,varargin);

% Call optimization function
switch class
    case 'skp'
        for i = 1:length(p)
            for j = 1:length(k)
                for l = 1:length(s)
                    [R,gamma] = Rskp(s,k,p);
                end
            end
        end
    
    case 'kp_imp'
        for i = 1:length(p)
            for j = 1:length(k)
                [R,alpha,beta]=Rkp_imp(k,p);
            end
        end
    
    case 'kp_dw'
        for i = 1:length(p)
            for j = 1:length(k)
                [R,alpha,beta,tbeta]=Rkp_dw(k,p);
            end
        end
    
    case 'kp_imp_dw'
        for i = 1:length(p)
            for j = 1:length(k)
                [R,alpha,beta,tbeta]=Rkp_imp_dw(k,p);  
            end
        end
end


multi_R = 1;
end



% =========================================================================

function [s]= setup_params(p,optional_params)
%function [s]= setup_params(p,optional_params)
%
% Set param values

i_p = inputParser;
i_p.FunctionName = 'setup_check_params';

% Default values
default_s = [];
default_writeToFile = 1;
    
% Populate input parser object
% ----------------------------
% Parameter values
i_p.addParamValue('s',default_s,@isnumeric);
i_p.addParamValue('writeToFile',default_writeToFile,@ischar);

i_p.parse(optional_params{:});
    
s = i_p.Results.s;    
end
% =========================================================================