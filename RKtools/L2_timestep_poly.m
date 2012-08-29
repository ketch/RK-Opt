function c = L2_timestep_poly(sdisc,p,q,doplot,eps,tol)
% function c = L2_timestep_poly(sdisc,p,q,eps,tol)
%
% Find the absolutely timestep for a given combination of
% linear spatial discretization and stability function.
%
% Also (optionally) plots the region of absolute stability and the eigenvalues.
%
% The timestep is determined to within accuracy eps (default 10^-4).
%
% The spectral stability condition is checked to within tol (default 10^-13).

if nargin<6  tol=1.e-13; end
if nargin<5  eps=1.e-5; end
if nargin<4  doplot=1; end
if nargin<3  q=zeros(size(p)); q(1)=1; end

if isfield(sdisc,'nx')==0 sdisc.nx=10; end

% Get semi-discrete operator
L = semispectrum(sdisc.name,sdisc.order,0,sdisc.nx);
lambda=eig(L);

prev=p(end:-1:1);
mv=0; cmin=0.; cmax=20.;
while cmax-cmin>eps
    c=(cmax+cmin)/2.;
    %Hacked for the moment:
    pval=ones(size(lambda));
    for ii=2:length(p)
        pval=pval+p(ii)*(c*lambda).^(ii-1);
    end
    mv=max(abs(pval));
    if mv>1+tol
        cmax=c;
    else
        cmin=c;
    end
    %Uncomment to unhack:
    %mv=max(abs(polyval(prev,c*lambda)./polyval(wrev(q),c*lambda)));
end
c=cmin;

if doplot
    %  figure
    %Determine region to plot
    plotbounds(1)=min(real(c*lambda))-4; plotbounds(2)=max(real(c*lambda))+1;
    plotbounds(3)=min(imag(c*lambda))-1; plotbounds(4)=max(imag(c*lambda))+1;
    %Plot absolute stability region and spectrum
    plotstabreg_func(p,q,plotbounds); hold on;
    semispectrum(sdisc.name,sdisc.order,doplot,sdisc.nx,c);
    
    %Plot roots
    zer=roots(prev);
    plot(real(zer),imag(zer),'ok','markersize',10)
    
    hold off; drawnow;
end
