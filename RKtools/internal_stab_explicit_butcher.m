function [stability] = internal_stab_explicit_butcher(A,b,c,spectrum,one_step_dt)
%function [stability] = internal_stability(rk_coeffs_file,rk_form,spectrum_file)
%
%
% This function compute and plots both intermediate and internal stability 
% functions of an explicit Runge-Kutta scheme given its Butcher tableau.
%
% Note that for an explicit Runge-Kutta scheme the stability functions are
% polynomial in the complex variable z.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check intermediate stability functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number os stages
s = length(A);

% Construct the intermediate stability functions \psi_j (where j is the 
% index of the stage).
%
% Note that for an explicit scheme the intermediate stability polynomial 
% associated to the first stage is always 1, i.e. \psi_1 = 1.
% Therefore we just compute and plots the remaining (s-1) intermediate
% stability polynomials plus the one-step stability polynomial of the
% Runge-Kuatta method.

% figure name
intermediate_stablity = figure;

% Grid 
hand = gca;
x = linspace(-40,2,200);
y = linspace(-40,20,200);
[X,Y] = meshgrid(x,y);
Z = (X + Y*sqrt(-1));
ac = 'ciao'

cmap = hsv(s);
for i = 1:s-1
    % Get generalized A and b^T
    [gA,gbT] = sub(A,i)
    
    % Cosntruct identity matrix of the right size
    Id = eye(length(gA));
    
    % Construct unit vector of the right size
    ge = ones(length(gA),1);
    
    % Define symbolic stability function
    syms z
    stab_function = matlabFunction(1 + gbT*z*(Id - gA*z)^(-1)*ge);
    
    [delta_t_inter] = stable_time_step(stab_function,1);
    
    % Evaluate stability function and take its module
    stab_function_eval = abs(subs(stab_function,{z},{Z}));
    
    %sprintf('%3u',i+1)
    i_stage = num2str(i+1);
    name = strcat(i_stage,' stage');
    % Plot contour   
    contour(X,Y,stab_function_eval,[1 1],'Color',cmap(i,:),'LineWidth',2,'DisplayName',sprintf(name));
    axis auto;
    hold on;
    
end

% Now plot the classical one-step stability function
Id =  eye(length(A));
e = ones(length(A),1);
syms z
stab_function = matlabFunction(1 + b'*z*(Id - A*z)^(-1)*e);
stab_function_eval = abs(subs(stab_function,{z},{Z}));
contour(X,Y,stab_function_eval,[1 1],'Color',cmap(s,:),'LineWidth',2,'DisplayName','One-step stability polynomial');
lg = legend('show');
set(lg,'Location','Best');
hold off;



end


% =========================================================================


% =========================================================================


function [delta_t] = stable_time_step(stab_fun,spec)
%function 
%
% Bisection method to find the stable time step for a given stability
% function
%%%%%%%%%%%%%%%%%%%%%%%%

% Bisection parameters
eps = 1.e-6;
toll = 1.e-10;

delta_t_min = 0.0; 
delta_t_max = 20.0;
max_z = 0.0;

% Symbolic independent variable
syms z

while (delta_t_max - delta_t_min) > eps
    delta_t = (delta_t_min + delta_t_max)/2.0;
    
    scaled_spec = delta_t*spec;
    
    stab_fun_eval = subs(stab_fun,{z},{scaled_spec});
    
    max_Z = max(abs(stab_fun_eval));
    
    if max_z > (1 + toll)
        delta_t_max = delta_t;
    else
        delta_t_min = delta_t;
    end
end

end
% =========================================================================


% =========================================================================

function [C,w] = sub(K,i)
%function mat = sub(K,i)
%
% Extract from the matrix K, the leading principal submatrix of order i
% K_i = K(1:i,1:i), and the vector of the weights.
%%%%%%%%%%%%%%%%%%%%%%%%

C = K(1:i, 1:i);
w = K(i+1,1:i);

% =========================================================================

end
