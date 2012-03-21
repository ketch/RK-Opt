function r = am_radius(rk,eps,rmax)
%begin_html
%function r = am_radius(rk,eps,rmax)
%
%By David Ketcheson
%
%Evaluates the Radius of absolute monotonicity
%of a Runge-Kutta method, given the Butcher array.
%
%For an m-stage method, rk.A should be an m x m matrix
%and rk.b should be a column vector of length m.
%
%Accuracy can be changed by modifying the value of eps (default 10^-10)
%Methods with very large radii of a.m. (>50) will require
%the default value of rmax to be increased.
%
% The radius of absolute monotonicity is the largest value of $r$
% such that
%\begin{eqnarray}
%K(I+rA)^{-1} &     \ge & 0 \\
%rK(I+rA)^{-1}e_m & \le & e_{m+1} 
%\end{eqnarray}
%where $$ K = \left(\begin{array}{c} A \\ b^T \end{array}\right) $$
%end_html


if nargin<3 rmax=50; end
if nargin<2 eps=1.e-10; end

m=length(rk.b); e=ones(m,1);
K=[rk.A;rk.b'];
rlo=0; rhi=rmax;

while rhi-rlo>eps  %use bisection
  r=0.5*(rhi+rlo);
  X=eye(m)+r*rk.A; beta=K/X; ech=r*K*(X\e);
  if (min(beta(:))<-3.e-16 || max(ech(:))>1.+3.e-16)
    rhi=r;
  else
    rlo=r;
  end
end

if rhi==rmax % r>=rmax
  error('Error: increase value of rmax in am_radius.m');
else
  r=rlo;
end
