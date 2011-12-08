********************************************************************************
Optimization of RK methods from multi-stability functions
********************************************************************************

The MATLAB scripts and functions contained in this directory extend
the original RK-opt/general package which is designed to search for optimized 
Runge-Kutta methods. These scripts optimize the Runge-Kutta coefficients 
starting from single or multi-stability functions. They also enable to use 
easily the MATLAB global optimization toolbox on multicore machines. The 
multistart solver is used to return local and global minima using different
starting points.

To run in parallel the user just need to run the command  "matlabpool open N" 
which opens N pool of MATLAB sessions for parallel computation, and execute
grk_opt_1stabPoly.m (single stability function) or grk_opt_mStabPoly.m
(multi-stability functions) with the correct inputs.


