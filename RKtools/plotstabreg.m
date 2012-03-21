function plotstabreg(rk,plotbounds,ls,lw)
%begin_html
%function plotstabreg(rk,plotbounds,ls,lw)
%
%Author: David Ketcheson
%
%Plots the absolute stability region
%of a Runge-Kutta method, given the Butcher array

%Inputs:
%  rk: a Runge-Kutta method
%Remaining inputs are optional:
%  plotbounds: bounds for region to compute and plot (default [-9 1 -5 5])
%  ls:   line style (default '-r')
%  lw:   line width (default 2)
%end_html
if nargin<4 lw=2; end
if nargin<3 ls='-r'; end
if nargin<2 plotbounds=[-9 1 -5 5]; end

%begin_html 
%Compute the stability function $\phi$ 
%end_html
[p,q]=rk_stabfun(rk);

plotstabreg_func(p,q,plotbounds,ls,lw);
