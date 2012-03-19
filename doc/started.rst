.. _installation:


===============
Installation
===============
This page describes how to obtain and test RK-opt. 

Dependencies 
------------
RK-opt depends on MATLAB 7.X or greater, and also requires the MATLAB optimization
toolbox. In order to 
test the basic functinalities of the code, we recommend that you also install the 
MATLAB xUnit test framework: `<http://www.mathworks.com/matlabcentral/fileexchange/22846>`_
MATLAB xUnit can be used with MATLAB R2008a or later. MATLAB xUnit relies 
heavily on object-oriented language features introduced in R2008a and will not 
work with earlier releases.

You will also need a git client on your system to obtain RK-Opt. itself


MATLAB xUnit test framework
+++++++++++++++++++++++++++
xUnit test framework can be downloaded for free at `<http://www.mathworks.com/matlabcentral/fileexchange/22846>`_
(click on the upper right corner button). 
An easy way to install xUnit without setting any environment variables is
to set up a MATLAB **startup.m** file, which executes commands of 
your choosing when the MATLAB program starts. 

If you want to follow this approach, you should add to your startup.m file the
following lines:

.. doctest::

   >>>> addpath "path-to-matlab_xunit_3"/matlab_xunit
   >>>> addpath "path-to-matlab_xunit_3"/matlab_xunit/xunit
   >>>> addpath "path-to-RK-opt"/RK-coeff-opt


Obtaining RK-opt
------------------
A zip file containing the package can be downloaded from
https://github.com/ketch/RK-opt/downloads.

If you wish to contribute back (and we hope you will), we recommend that you
fork the RK-Opt Github repository, implement your additions, and issue
a pull request.  You may also simply e-mail a patch to us.


Testing your installation with MATLAB xUnit test framework
----------------------------------------------------------
You can test your RK-opt installation with xUnit.  First ::

    $ cd "path-to-RK-Opt"/general/test

Then run the following command in the GUI

.. doctest::

    >>> runtests

If everything is set up correctly, this will run several tests, and inform you 
that the tests passed. Although the tests available are not very extensive, they
cover some of the main functions of the package. An extension of the test suite
is planned.
