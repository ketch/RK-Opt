.. _codeStructure:


=================
RK-Opt. structure
=================
The main function of the code is **rk_opt** which is stored in **rk_opt.m**. 
Given the coefficients of a stability function :math:`R(z)`, the order of 
accuracy :math:`p` of the scheme and its number of stages :math:`s`, 
rk_opt setups and calls the MATLAB's fmincon routine to find a Runge-Kutta 
method that satisfies both equality constrains and a combination of inequality 
constraints and objective functions.


rk_opt function signature 
-------------------------
The **rk_opt** function signature is ::

     rk_opt(s,p,class,objective,poly_coeff_ind,poly_coeff_val,startvec,solveorderconditions,np,max_tries,writeToFile).

The meaning of the arguments is as follow:
    * :math:`s`: number of stages
    * :math:`p`: order of the Runge-Kutta (RK) scheme
    * class: class of method to search (erk = explicit RK; irk = implicit RK; dirk = diagonally implicit RK; sdirk = ingly diagonally implicit RK; 2S, 3S, 2S*, 3S* = low-storage formulations, see *Ketcheson, "Runge-Kutta methods with minimum storage implementations". J. Comput. Phys. 229(5):1763 - 1773, 2010*)
    * objective: objective function (acc = minimize leading truncation error coefficient; ssp = maximize SSP coefficient)
    * poly_coeff_ind: index of the polynomial coefficients for :math:`j > p`  
    * poly_coeff_val: of the polynomial coefficients for :math:`j > p` (tall-tree elementary weights)
    * startvec: vector of the initial guess
    * solveorderconditions: if set to 1, solve the order conditions first before trying to optimize
    * np: number of processor to use. If :math:`np > 1` the MATLAB global optimization toolbox *Multistart* is used
    * max_tries: maximum number of **fmincon** function calls
    * writeToFile: wether to write to a file. If set to 1 write the RK coefficients to a file called "ERK-p-s.txt"

    


   
