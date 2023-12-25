% Code to look for PEP methods
% stages and classical order of accuracy
s = 6
p = 2
% Number of cores (4 in my computer)
n_cores = 4

% This line was recommended by David: it creates a loop which looks for methods satisfying the nonlinear conditions.
% We want to make 1000 searches with each algorithm (sqp, interior-point, active-set). 
% If no method is found, then we can move to the next number of stages.
rk = rk_opt(s,p,'erk','acc', 'algorithm', 'sqp','np',6,'solveorderconditions',1,'num_starting_points',1500,'suppress_warnings',1);

% I think the the last command may run forever until it finds a solution (correct me if I am wrong). Since I don't 
% know how to keep track of the number of searches, I made another loop. I can include it later in this file
% if necessary. 


