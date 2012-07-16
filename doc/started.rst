.. _installation:


===============
Installation
===============
This page describes how to obtain and test RK-opt. 

Dependencies 
------------
 - MATLAB 7.X or greater
 - MATLAB optimization toolbox

Optional dependencies
------------------------------
 - MATLAB Global Optimization toolbox (for multicore search)
 - xUnit test framework (for testing your installation); available for free
   from `<http://www.mathworks.com/matlabcentral/fileexchange/22846>`_.
   xUnit requires MATLAB R2008a or later. 


Obtaining RK-opt
------------------
 - Download: https://github.com/ketch/RK-opt/downloads.  
 - Or clone::

    $ git clone https://github.com/ketch/RK-opt.git

After unzipping/cloning, add the subdirectory `RK-opt/RKtools` to your MATLAB path.


=========================
Testing your installation 
=========================
You can test your RK-opt installation with xUnit.  

Installing xUnit
----------------
The MATLAB xUnit test framework can be downloaded for free at
`<http://www.mathworks.com/matlabcentral/fileexchange/22846>`_
(after the link, click on the button in the upper right corner). 
An easy way to install xUnit without setting any environment variables is
to add the following lines to your **startup.m** file:

.. doctest::

   addpath "path-to-matlab_xunit_3"/matlab_xunit
   addpath "path-to-matlab_xunit_3"/matlab_xunit/xunit
   addpath "path-to-RK-opt"/RK-coeff-opt

Running the tests
-----------------

To run the tests, do the following in MATLAB::

    >>> cd "path-to-RK-Opt"/general/test
    >>> runtests

If everything is set up correctly, this will run several tests, and inform you 
that the tests passed. At present the tests are not very extensive.
