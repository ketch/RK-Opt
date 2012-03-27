function multi_R = multi_R_opt(k,p,class,varargin)
%multi_R = multi_R_opt(k,p,class,varargin)
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
% when optimal contractive k-step, s-stage GLM are investigated.
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
[s]= setup_params(p,class,varargin);

% Set output file name
output_file_name = strcat(class,'.txt');

% Open output file
write_fid = fopen(output_file_name,'w');

% Write to file number of optimization runs
number_runs = write_number_runs(k,p,s,class,write_fid);


% Call optimization function
switch class
    case 'skp'
        for i = 1:length(p)
            for j = 1:length(k)
                for l = 1:length(s)
                    glm.order = p(i);
                    glm.step = k(j);
                    glm.stage = s(l);
                    [glm.R,glm.gamma] = Rskp(glm.stage,glm.step,glm.order);
                    glm.coeffs_monomial = gamma2monomial_coeffs(glm.gamma,glm.R);
                    wf = write_output_to_file(write_fid,glm);
                end
            end
        end
    
    case 'kp_imp'
        for i = 1:length(p)
            for j = 1:length(k)
                glm.order = p(i);
                glm.step = k(j);
                [glm.R,glm.alpha,glm.beta]=Rkp_imp(glm.step,glm.order);
                wf = write_output_to_file(write_fid,glm);
            end
        end
    
    case 'kp_dw'
        for i = 1:length(p)
            for j = 1:length(k)
                glm.order = p(i);
                glm.step = k(j);
                [glm.R,glm.alpha,glm.beta,glm.tbeta]=Rkp_dw(glm.step,glm.order);
                wf = write_output_to_file(write_fid,glm);
            end
        end
    
    case 'kp_imp_dw'
        for i = 1:length(p)
            for j = 1:length(k)
                glm.order = p(i);
                glm.step = k(j);
                [glm.R,glm.alpha,glm.beta,glm.tbeta]=Rkp_imp_dw(glm.step,glm.order); 
                wf = write_output_to_file(write_fid,glm);
            end
        end
end

multi_R = 1;

end

% =========================================================================


% =========================================================================

function [s]= setup_params(p,class,optional_params)
%function [s]= setup_params(p,class,optional_params)
%

% Set param values
i_p = inputParser;
i_p.FunctionName = 'setup_check_params';

% Default values
default_write_file = 1;
    
% Populate input parser object
% ----------------------------
% Parameter values
if strcmp(class,'skp')
    default_s = p;
else
    default_s = [];
end

i_p.addParamValue('s',default_s,@isnumeric);
i_p.addParamValue('write_file',default_write_file,@ischar);

i_p.parse(optional_params{:});
    
s = i_p.Results.s;    
end
% =========================================================================


% =========================================================================

function number_runs = write_number_runs(k,p,s,class,write_fid)
%function  number_runs = write_number_runs(k,p,s,class,write_fid)
%
if strcmp(class,'skp')
    number_runs = length(p)*length(k)*length(s);
else
    number_runs = length(p)*length(k);
end

str = 'stability polynomials';
fprintf(write_fid,'%s\r\n',str);

fprintf(write_fid,'%u\n',number_runs);
str = '=======================================================================';
fprintf(write_fid,'%s\r\n',str);

end
% =========================================================================


% =========================================================================

function wf = write_output_to_file(write_fid,glm)
%function  wf = write_output_to_file(write_fid,glm)
%

% Convert glm to cell and extract values and names
values = struct2cell(glm);
names  = fieldnames(glm);

% Write field
for i=1:length(values)
    write_field(write_fid,names{i},values{i});
end

str = '=======================================================================';
fprintf(write_fid,'\n%s\r\n\n',str);


wf = 1;
end
% =========================================================================


% =========================================================================

function wf = write_field(write_fid,name,value)
%function wf=writeField(write_fid,name,value)
%

fprintf(write_fid,'\n%s\n',name);
fprintf(write_fid, [repmat('%5.16E\t', 1, size(value,1)),'\n'], value');

wf = 1;
end
% =========================================================================


% =========================================================================

function a = gamma2monomial_coeffs(gamma,R)
%function a = gamma2monomial_coeffs(gamma,R)
%
% gamma = coefficients of the stability polynomials using the basis (1+z/R)^j 
%
% a = coefficients in the monomial basis
% 
% Transformation: 
%    
%    a_i = 1/(i!*R^i) * sum_{j=0}^m (gamma_j prod_{n=0}^{i-1}(j-n))
%

for i = 0:length(gamma)-1 % loop over coefficients of the monomial basis
    sm = 0;
    for j = 0:length(gamma)-1 % summation over gamma_j
        pr = 1;
        for n = 0:i-1 % product (j-n)
            pr = pr*(j-n);
        end
        % The algorithm is limited by numerical conditioning of the 
        % constraint matrix. Nagative (hopefully very small) gamma_j are 
        % due to this limitation. We set them to zero.
        if (gamma(j+1) < 0.0)
            gamma(j+1) = 0.0;
        end
        sm = sm + gamma(j+1)*pr; 
    end
    a(i+1) = 1/(factorial(i)*R^i)*sm;
end

a = a'; % Transpose for pretty output

end
% =========================================================================



