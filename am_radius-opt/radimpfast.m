%van de Griend's algorithm
%Assuming multiplicity one for all roots
%Uses high precision arithmetic
function rad=radimpfast(p,q)

syms x P Q phi

digits(100);
maple('Digits:=100');
meps=1.e-8;
eps=1.e-10; M=20;

rzero=false; m=length(p)-1; n=length(q)-1;

% evaluate phi = P/Q
i=1:length(p); P=sum(p(i).*x.^(i-1));
i=1:length(q); Q=sum(q(i).*x.^(i-1));
phi=P/Q

% Compute the partial fractions decomposition of phi
phi=maple('convert',phi,'parfrac','complex');
phi2=maple('convert',phi,'list');
sd=0;
for i=1:m-n+1 %Peel off direct terms
  sd(i)=maple('coeff',phi2(i),'x',m-n-i+1);
end
tt=max(m-n+1,0); ii=1;
for i=tt+1:length(phi2) %peel off poles
  spoles(ii)=-maple('op',2,maple('denom',phi2(i)))/maple('coeff',maple('denom',phi2(i)),'x',1);
  sc(ii)=-maple('numer',phi2(i)/maple('coeff',maple('denom',phi2(i)),'x',1));
  ii=ii+1;
end
d=double(sd);dpoles=double(spoles);c=double(sc);
lk=length(d);
npoles=length(dpoles);

[alpha0,ind]=min(abs(dpoles));
if abs(imag(dpoles(ind)))>meps || real(dpoles(ind))<=0
  rzero=true;
end
%Reorder so alpha0 is first
rpoles=dpoles; rspoles=spoles;
rpoles(1)=dpoles(ind); rspoles(1)=spoles(ind);
rpoles(ind)=dpoles(1); rspoles(ind)=spoles(1);
rc=c; rsc=sc;
rc(1)=c(ind); rsc(1)=sc(ind);
rc(ind)=c(1); rsc(ind)=sc(1);

%Find B
B=M;
for i=2:npoles
  if (abs(imag(rpoles(i)))>meps || real(rpoles(i))<0)
    B=min(B,-double((abs(rspoles(i))^2-alpha0^2)/(2*(real(rspoles(i))-alpha0))));
  end
end
if B<0 rzero=true; end
if rc(1)<0 rzero=true; end

%Find K(0)
if rzero==false

  L0=max((alpha0-abs(rpoles(2:npoles)))./(abs(rpoles(2:npoles))-alpha0));
  L0=max([L0,0,m-n+1]);
  k=ceil(L0-1); F0=rc(1)*10;
  while rc(1)<F0
    k=k+1;
    F0=sum(abs(rc(2:npoles)).*(abs(alpha0)^(k+1)./abs(rpoles(2:npoles)).^(k+1)));
  end
  K0=k;k=0;phid=phi;
  qhat=(-1).^(0:length(q)-1).*q;
  while (k<K0 && rzero==false)
    if k>0 phid=simplify(diff(phi)); end
    if k<=max(m,n)
      val=double(maple('eval',phid,'x=0'));
      beta(k+1)=(-1)^k*val/factorial(k);
    else
      beta(k+1)=-1./qhat(1)*sum(qhat(2:n+1).*beta(k:-1:(k-n+1)));
%      beta(k+1)*(-1)^k*factorial(k)-subs(diff(phi,k),0)
    end
    if beta(k+1)*(-1)^k<=0 rzero=true; end
    k=k+1;
  end
end

if rzero==false
  rlo=max(B-eps,0);
  x0=-rlo;
  L=max((abs(alpha0-x0)-abs(rpoles(2:npoles)-x0))./(abs(rpoles(2:npoles)-x0)-abs(alpha0-x0)));
  L=max([L,0,m-n+1]); K=ceil(L-1); F=rc(1)*10;
  while rc(1)<F
    K=K+1;
    F=sum(abs(rc(2:npoles)).*(abs(alpha0-x0)^(K+1)./abs(rpoles(2:npoles)-x0).^(K+1)));
  end
  kch=0; phid=phi;

  while kch<K+1 %Check that first K derivatives are positive
%    if kch>0 phid=simplify(diff(phid)); end %Get next derivative
%    if kch>0 phid=maple('simplify',(maple('diff',phid,x))); end %Equally fast
%    if kch>0 phid=maple('diff',phi,['x$',num2str(kch)]); end
%      val=real(double(maple('eval',phid,['x=',num2str(x0,14)])))
      if lk>kch %still have direct terms
        val=0;
        for ii=1:lk
          val=val+sd(ii).*x0.^(lk-ii-kch).*prod(lk-ii+1-kch:lk-ii);
	end
	val=val+(-1)^(kch+1)*factorial(kch)*sum(sc./(x0-spoles).^(kch+1));
      else
        val=    (-1)^(kch+1)*factorial(kch)*sum(sc./(x0-spoles).^(kch+1));
      end
      val=real(double(val));
%      val=real(subs(phid,x0));
      if val<-1.e-14 %Phi is not a.m. here
        xhi=-x0; xlo=0;
        while xhi-xlo>eps
          x1=-0.5*(xhi+xlo);
          if lk>kch %still have direct terms
            v=0;
            for ii=1:lk
              v=v+sd(ii).*x1.^(lk-ii-kch).*prod(lk-ii+1-kch:lk-ii);
            end
	    v=v+(-1)^(kch+1)*factorial(kch)*sum(sc./(x1-spoles).^(kch+1));
          else
            v=  (-1)^(kch+1)*factorial(kch)*sum(sc./(x1-spoles).^(kch+1));
          end
	  v=real(double(v));
	  if (v<0 || isnan(v)) xhi=-x1; else xlo=-x1; end
        end %while
        x1=-xlo;
	if x1>x0 x0=x1; kch=-1; phid=phi;
        L=max((abs(alpha0-x0)-abs(rpoles(2:npoles)-x0))./(abs(rpoles(2:npoles)-x0)-abs(alpha0-x0)));
        L=max([L,0,m-n+1]);
        k=ceil(L-1); F=rc(1)*10;
        while rc(1)<F
          k=k+1;
          F=sum(abs(rc(2:npoles)).*(abs(alpha0-x0)^(k+1)./abs(rpoles(2:npoles)-x0).^(k+1)));
        end %while
        K=k;
        end %if x1>x0
      end %if val<0
    kch=kch+1;
  end
  rad=-x0;
%  rad=B;
else 
  rad=0;
end
rad=-rad;
