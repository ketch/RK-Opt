function multi_rk = multi_am_radius_2_rk_opt(inputFileName,class,varargin)
%function multi_rk = multi_am_radius_2_rk_opt(inputFileName,class,varargin)
%
%
% This function calls the rk_opt.m function inside a loop for optimizing
% multiple explicit monotonicity preserving RK methods given their stability 
% function coefficients. 
%
% Currently the stability function coefficients in rk_opt.m is imposed using the
% monomial basis, i.e. \Psi(z) = sum_{j=0}^m a_j z^j.
%
%
% inputFileName: name of the file where the coefficients of the stability 
%                polynomial function can be found. This file must have a
%                specific structure. For instance:
%   
%                #stability poly.
%                 1
%
%                 ==========================================================
%
%
%                 order
%                 2.0000000000000000E+00	
%
%                 step
%                 1.0000000000000000E+00	
%
%                 stage
%                 2.0000000000000000E+00	
%
%                 R
%                 1.0000000099992548E+00	
%
%                 gamma
%                 5.0000000000000011E-01	-9.9992548323769717E-09	5.0000000999925476E-01	
%
%                 coeffs_monomial
%                 1.0000000000000000E+00	1.0000000000000000E+00	5.0000000000000000E-01	
%
%                 ==========================================================
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% =========================================================================
% Load and read the file containing the stability polynomial coefficients
% =========================================================================
