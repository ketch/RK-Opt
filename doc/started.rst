.. _installation:


===============
Installation
===============
This page describe how to install and test RK-opt. The installation is extremely 
simple. However, if you find any diffuculties during the process, please send an 
e-mail to xxxxxxxxx.

Dependencies 
------------
RK-opt depends for the moment on MATLAB 7.X or greater. However, in order to 
test the basic functinalities of the code, we recommend to also install the 
MATLAB xUnit test framework: `<http://www.mathworks.com/matlabcentral/fileexchange/22846>`_
MATLAB xUnit can be used with MATLAB R2008a or later. MATLAB xUnit relies 
heavily on object-oriented language features introduced in R2008a and will not 
work with earlier releases.

You will also need a git client on your system to obtain RK-Opt. itself.


Obtaining MATLAB?
+++++++++++++++++
MATLAB is not free and the act of taking and using software without having paid
for it is piracy and a crime. Since we would like to give you free 
access to our code without violate any law, we are planning to have a 
Python-based version of RK-opt. 

MATLAB xUnit test framework
+++++++++++++++++++++++++++
xUnit test framework can be downloaded for free at `<http://www.mathworks.com/matlabcentral/fileexchange/22846>`_
(click on the upper right corner button). 
An easy way to install xUnit without setting up any variable in your shell 
profile is to setup a MTALB **startup** file. startup executes commands of 
your choosing when the MATLAB program starts. 

If you want to follow this approach, you should add to your startup.m file the
following lines:

.. doctest::

   >>>> addpath "path-to-matlab_xunit_3"/matlab_xunit
   >>>> addpath "path-to-matlab_xunit_3"/matlab_xunit/xunit


Installing RK-opt
------------------
The current method for installing RK-opt is to create a local copy of its 
github-hosted repository::

    $ git@github.com:ketch/RK-opt.git


Testing your installation with MATLAB xUnit test framework
----------------------------------------------------------

You can test your RK-opt installation with xUnit. Thus ::

    $ cd "path-to-RK-Opt"/general/test

Then run the following command in the GUI

.. doctest::

    >>> runtests

If everything is set up correctly, this will run several tests, and inform you 
that the tests passed. Although the tests available are not very extensive, they
cover some of the main functions of the package. An extension of the tests suite
is planned in the near future.



