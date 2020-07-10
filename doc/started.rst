.. _installation:


===============
Installation
===============
This section describes how to obtain RK-Opt and test that it is working correctly.

Dependencies
------------
 - MATLAB (relatively recent versions; tested with R2018a and later)
 - MATLAB Optimization Toolbox

Optional dependencies
------------------------------
 - MATLAB Global Optimization toolbox (for multicore search)


Obtaining RK-Opt
------------------
 - Download: https://github.com/ketch/RK-Opt/
 - Or clone::

    $ git clone https://github.com/ketch/RK-Opt.git

After unzipping/cloning, add the subdirectory ``RK-Opt/RKtools`` to your MATLAB path.


=========================
Testing your installation
=========================
You can test your RK-Opt installation by running the MATLAB script `test.m`.

Running the tests
-----------------

To run the tests, do the following in MATLAB::

    >> cd /path/to/RK-Opt/
    >> test

If everything is set up correctly, this will run several tests, and inform you
that the tests passed.
