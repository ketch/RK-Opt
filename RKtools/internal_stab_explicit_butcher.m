function [stability] = internal_stab_explicit_butcher(A,b,c,p,spectrum,one_step_dt)
% function [stability] = internal_stab_explicit_butcher(A,b,c,spectrum,one_step_dt,p)
%
%
% This function computes and plots both intermediate and one-step internal 
% stability vector of an explicit Runge-Kutta scheme given its Butcher 
% tableau.
%
% Note that for an explicit Runge-Kutta scheme the stability functions are
% polynomials in the complex variable z.
%
% Construct the intermediate stability functions \psi_j (where j is the 
% index of the stage).
%
% Note that for an explicit scheme the intermediate stability polynomial 
% associated to the first stage is always 1, i.e. \psi_1 = 1.
% Therefore we just compute and plot the remaining (s-1) intermediate
% stability polynomials plus the one-step stability polynomial of the
% Runge-Kuatta method.


% Number of stages
s = length(A);

close all;
clc;

% Figure name
intermediate_stability = figure;

% Grid 
hand = gca;
x = linspace(-40,2,200);
y = linspace(-40,20,200);
[X,Y] = meshgrid(x,y);
Z = (X + Y*sqrt(-1));

% Create map of colors
cmap = hsv(s);
for j = 1:s-1
    % Get generalized A and b^T
    [gA,gbT] = sub(A,j);
    
    % Cosntruct identity matrix of the right size
    Id = eye(length(gA));
    
    % Construct unit vector of the right size
    ge = ones(length(gA),1);
    
    % Define symbolic stability function
    z = sym('z');
    stab_function = 1 + gbT*z*(Id - gA*z)^(-1)*ge;
    
    [delta_t_inter] = stable_time_step(stab_function,1);
    
    % Evaluate intermediate stability function and compute its module
    stab_function_eval = abs(subs(stab_function,{z},{Z}));
    
    % Plot contour   
    stage = num2str(j+1);
    name = strcat(stage, ' stage');
    contour(X,Y,stab_function_eval,[1 1],'Color',cmap(j,:),'LineWidth',2,'DisplayName',sprintf(name));
    axis auto;
    hold on;
    
end

% Construct and plot the classical one-step Runge-Kutta stability function
Id =  eye(length(A));
e = ones(length(A),1);
z = sym('z');
stab_function = 1 + b'*z*(Id - A*z)^(-1)*e;
stab_function_eval = abs(subs(stab_function,{z},{Z}));
contour(X,Y,stab_function_eval,[1 1],'Color',cmap(s,:),'LineWidth',2,'DisplayName','One-step stability polynomial');
lg = legend('show');
set(lg,'Location','NorthWest','FontSize',18,'FontWeight','bold');

xlabel('Re','FontSize',18,'FontWeight' ,'bold');
ylabel('Im','FontSize',18,'FontWeight' ,'bold');
set(hand,'FontSize',18,'FontWeight' ,'bold');

title('Intermediate stability polynomial','FontSize',18,'FontWeight' ,'bold')

grid on;

hold off;




% s-stage INTERNAL STABILITY VECTOR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute scaled spectrum
scaled_spectrum = spectrum*one_step_dt;

% Compute the s-stage internal stability vector (\theta).
Id =  eye(length(A));
z = sym('z');
theta = b'*z*(Id - A*z)^(-1);

% Evaluate the internal stability vector. Note that the first component of
% the vector is not used because for explicit scheme in first stage no
% round-off error is introduced.
theta_eval = abs(subs(theta(2:s),{z},{scaled_spectrum}));

% Compute the maximum norm
n_s= max(theta_eval);
wf = write_file(s,p,n_s);

% Plot contour of the s-stage stability vector components 
internal_stability = figure;
hand = gca;
for j = 2:s
    theta_stage = theta(j);
    theta_eval = abs(subs(theta_stage,{z},{Z}));
    stage = num2str(j);
    name = strcat(stage, ' stage');
    contour(X,Y,theta_eval,[1 1],'Color',cmap(j-1,:),'LineWidth',2,'DisplayName',sprintf(name));
    axis auto;
    hold on;
end
lg = legend('show');
set(lg,'Location','NorthWest','FontSize',18,'FontWeight','bold');

xlabel('Re','FontSize',18,'FontWeight' ,'bold');
ylabel('Im','FontSize',18,'FontWeight' ,'bold');
set(hand,'FontSize',18,'FontWeight' ,'bold');

title('One-step internal stability polynomials','FontSize',18,'FontWeight' ,'bold')

grid on;

hold off;

end
% =========================================================================


% =========================================================================
function [C,w] = sub(K,i)
%function [C,w] = sub(K,i)
%
% Extract from the matrix K, the leading principal submatrix of order i
% K_i = K(1:i,1:i), and the vector of the weights.
%%%%%%%%%%%%%%%%%%%%%%%%

C = K(1:i, 1:i);
w = K(i+1,1:i);

end
% =========================================================================


% =========================================================================
function [delta_t] = stable_time_step(stab_fun,spec)
%function [delta_t] = stable_time_step(stab_fun,spec)
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
function wf = write_file(s,p,val)
%function wf = write_file(s,p,n)
%
% 
% Write to file the maximum norm of the j-stage internal stability vector
% \theta_j, 1=< j <= s

output_file_name = strcat('max-norm_int-stab-',num2str(p),num2str(s),...
                          '.txt');
write_fid = fopen(output_file_name,'w');

fprintf(write_fid, '%s%u\n%s%u\n%s\t%3.7f\n','stage ',s,'order ',p,...
        'maximum norm ',val);

wf = 1;
end
% =========================================================================
