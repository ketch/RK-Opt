.. _RK-opt_tutorial:

********
Tutorial
********
As an example, in this tutorial we will use the RK-opt package to compute the
Butcher's tableau given the stability function of a 4th-order RK scheme optimized
for a 4th-order Spectral Difference spatial discretization.

PROVIDE THE TXT FILE. 


Single start
============
You can run the rk_opt function directly from the MATLAB GIU in this way ::

    >>>> rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,np,max_tries,writeToFile)


Multistart
==========
You can run the rk_opt function directly from the MATLAB GIU in this way ::

    >>>>multi_rk = multi_rk_opt(inputFileName,class,objective,startvec,solveorderconditions,np,max_tries,writeToFile)






