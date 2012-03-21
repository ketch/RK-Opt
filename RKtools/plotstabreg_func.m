function [dummy] = plotstabreg_func(p,q,bnds,ls,lw)
%begin_html
%function [dummy] = plotstabreg_func(p,q,bnds,ls,lw)
%
%plot the absolute stability region of a one-step method,
%given the stability function
%
%Inputs:
%  p: coefficients of the numerator   of the stability function
%  q: coefficients of the denominator of the stability function 
%  if q is omitted, it is assumed that the function is a polynomial
%Remaining inputs are optional:
%  bnds: bounds for region to compute and plot (default [-9 1 -5 5])
%  ls:   line style (default '-r')
%  lw:   line width (default 2)
%end_html


if nargin<5 lw=2; end
if nargin<4 ls='-r'; end
if nargin<3 bnds=[-9 1 -5 5]; end
if nargin<2 q=[1]; end

dx=(bnds(2)-bnds(1))/500.;
xa=bnds(1):dx:bnds(2); ya=bnds(3):dx:bnds(4);

X=repmat(xa,length(ya),1); Y=repmat(ya',1,length(xa));
XY=(X+1i*Y);
S=ones(size(XY)); T=ones(size(XY));
m=length(p)-1; n=length(q)-1;

%Evaluate numerator
XYP=XY;
for j=1:m
  S=S+p(j+1)*XYP; XYP=XYP.*XY;
end
%Evaluate denominator
XYP=XY;
for j=1:n
  T=T+q(j+1)*XYP; XYP=XYP.*XY;
end

%Compute $R(z)=|S(z)/T(z)|$
R=abs(S./T);

%Plot the absolute stability boundary ($R=1$)
contour(xa,ya,R,[1 1],ls,'LineWidth',lw);
title('Absolute stability region');
hold on
plot([0 0],bnds(3:4),'--k'); plot(bnds(1:2),[0 0],'--k')
axis(bnds)
hold off
