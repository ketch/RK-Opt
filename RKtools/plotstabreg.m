function [contour_matrix] = plotstabreg(rk,plotbounds,ls,lw)
% function [contour_matrix] = plotstabreg(rk,plotbounds,ls,lw)
%
% Plots the absolute stability region
% of a Runge-Kutta method, given the Butcher array

% Inputs:
%       * rk: a Runge-Kutta method

% Remaining inputs are optional:
%       * plotbounds: bounds for region to compute and plot (default [-9 1 -5 5])
%       * ls:   line style (default '-r')
%       * lw:   line width (default 2)

if nargin<4 lw=2; end
if nargin<3 ls='-r'; end
if nargin<2 plotbounds=[-9 1 -5 5]; end

[p,q]=rk_stabfun(rk);

contour_matrix = plotstabreg_func(p,q,plotbounds,ls,lw);

if isfield(rk, 'Ahat') && isfield(rk, 'bhat') && isfield(rk, 'chat')
    hold on
    rk_emb.A = rk.Ahat;
    rk_emb.b = rk.bhat;
    rk_emb.c = rk.chat;
    [phat, qhat] = rk_stabfun(rk_emb);
    plotstabreg_func(phat, qhat, plotbounds, '-b', lw);
    legend('Main', 'Embedded');
    hold off
end
