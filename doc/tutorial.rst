.. _RK-opt_tutorial:

********
Tutorial
********
As an example, in this tutorial we will use the RK-opt package to compute the
Butcher's tableau given the stability function of a 4th-order RK scheme optimized
for a 4th-order spatial discretization.
 

Single start
============
Data:
    * :math:`p` = 4
    * :math:`s` = 10
    * stability polynomial coefficients: 1.0000000000000000E+00	1.0000000000000000E+00	5.0000000000000000E-01	1.6666666666666666E-01	4.1666666666666664E-02	7.9381511986016205E-03	1.1433113866910084E-03	1.2103098081723528E-04	8.9101339095924001E-06	4.0846353688297018E-07	8.7981307960076615E-09

You can run the rk_opt function directly from the MATLAB GIU in this way ::

    >>>> rk = rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,np,max_tries,writeToFile)


Multistart
==========
Matteo will provide the .txt file with several stability polynomial coefficients 
sets.

You can run the rk_opt function directly from the MATLAB GIU in this way ::

    >>>> multi_rk = multi_rk_opt(inputFileName,class,objective,startvec,solveorderconditions,np,max_tries,writeToFile)






