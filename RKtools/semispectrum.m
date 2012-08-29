function L = semispectrum(method,order,doplot,nx,cfl)
% function L = semispectrum(method,order,doplot,nx,cfl)
% Plot spectra of various semi-discretizations of the advection equation
%
% Current choices for method:
%       - 'fourier':   Fourier   spectral method
%       - 'chebyshev': Chebyshev spectral method
%       - 'updiff':    Upwind difference operators (linearized WENO)
%       - 'DG':        Discontinuous Galerkin method
%
% The value of order matters only for the 'updiff' and 'DG' methods
% and selects the order of accuracy in those cases.

if nargin<5 cfl=1; end
if nargin<4 nx=100; end
if nargin<3 doplot=1; end

switch method

  case 'heat_cdiff' %Centered-difference for the heat equation
    L=-2*diag(ones(nx,1),0)+diag(ones(nx-1,1),-1)+diag(ones(nx-1,1),1); L(1,end)=1; L(end,1)=1;
    h=1./(nx+1);
    L=L/h^2;
  case 'fourier'
    column=[0 0.5*(-1).^(1:nx-1).*cot((1:nx-1)*pi/nx)];
    L=toeplitz(column,column([1 nx:-1:2]));
  case 'chebyshev'
    L=cheb(nx);
  case 'legendre'
    [x,L]=legDc(nx);
  case 'updiff'  %Linearized WENO (i.e. upwind difference) operators
    e=ones(nx,1);
    switch order
      case 1
        L=spdiags([-e,e],[0,-1],nx,nx); L(1,end)=1;
        L=full(L);
      case 3
        L2a=full(spdiags([e/2 e/2],[0 1],nx,nx));
        L2a(end,1)=L2a(end-1,end);
	L2b=full(spdiags([-e/2 e*3/2],[-1 0],nx,nx));
	L2b(1,end)=L2b(2,1);
	L2=2/3*L2a+1/3*L2b;
        L1=[L2(:,2:end) L2(:,1)];
        L=L1-L2;
      case 5
        L2a=full(spdiags([e/3 e*5/6 -e/6],[0 1 2],nx,nx));
        L2a(end,1)=L2a(end-1,end);
        L2a(end-1,1)=L2a(end-2,end); L2a(end,2)=L2a(end-2,end);
        L2b=full(spdiags([-e/6 e*5/6 e/3],[-1 0 1],nx,nx));
        L2b(1,end)=L2b(2,1);
        L2b(end,1)=L2b(end-1,end);
        L2c=full(spdiags([e/3 -e*7/6 e*11/6],[-2 -1 0],nx,nx));
        L2c(1,end)=L2c(2,1);
        L2c(2,end)=L2c(3,1); L2c(1,end-1)=L2c(3,1);
        L2=3/10*L2a+3/5*L2b+1/10*L2c;
        L1=[L2(:,2:end) L2(:,1)];
        L=L1-L2;
      otherwise
        k=(order+1)/2;
        coeffs=-fdcoeffF(1,k+1,1:2*k);
	M=e*coeffs;
	L=full(spdiags(M,-k:(k-1),nx,nx));
	for i=1:k
          L(i,(end-k+i):end)=coeffs(1:k-i+1);
	end
	for i=1:k-1
	  L(end-k+i+1,1:i)=coeffs(end-i+1:end);
	end
      end

  case 'DG'
    %Discontinuous Galerkin with periodic BCs
    Luw=-diag(ones(nx,1),0)+diag(ones(nx-1,1),-1); Luw(1,end)=1;
    L2=-diag(ones(nx,1),0)-diag(ones(nx-1,1),-1); L2(1,end)=-1;
    switch order
      case 1 % 1st order:
        L=Luw;
      case 2 % 2nd order:
        L=[Luw Luw;-3*Luw 3*L2];
      case 3 % 3rd order:
        L=[Luw Luw Luw;-3*Luw 3*L2 3*L2; 5*Luw -5*L2 5*Luw];
      case 4 % 4th order: (THIS IS WRONG!)
        L=[Luw Luw Luw Luw;-3*Luw 3*L2 3*L2 3*L2; 5*Luw -5*L2 5*Luw 5*Luw;
	-7*Luw 7*L2 -7*Luw 7*L2];
    end
  case 'compact'
    switch order
      case 10 %Tenth order compact scheme of AFLES: L*u^n = (Ra+nu*Rd)*u^{n-1}
        v=[1 1/2 1/20 zeros(1,nx-5) 1/20 1/2];
	L=toeplitz(v); %Left-hand-side
	a=17/12*ones(nx-1,1)/2.; b=101/150*ones(nx-2,1)/4.; c=1/100*ones(nx-3,1)/6.;
	aa=17/12*ones(1,1)/2.; bb=101/150*ones(2,1)/4.; cc=1/100*ones(3,1)/6.;
        Ra=diag(zeros(nx,1))+diag(a,1)-diag(a,-1)+diag(b,2)-diag(b,-2)...
            +diag(c,3)-diag(c,-3)...
            -diag(aa,nx-1)-diag(bb,nx-2)-diag(cc,nx-3)...
            +diag(aa,-nx+1)+diag(bb,-nx+2)+diag(cc,-nx+3);
	L=L\Ra;
    end

end

if doplot
  %Plot spectrum of semi-discrete scheme:
  z=eig(cfl*L); plot(real(z),imag(z),'.g','markersize',10);
  %axis equal
end
%ps_fun(L);
