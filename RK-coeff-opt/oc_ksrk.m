function coneq= oc_ksrk(A,b,D,theta,p)
% function coneq= oc_ksrk(A,b,D,theta,p)
% Order conditions for multistep-RK methods.

s=length(b);
k=length(theta);
es=ones(s,1);
ll=[k-1:-1:0]';
c=A*es-D*ll;
b=b';

% tau(k)= c.^k-D*(-ll').^k/k!- A*c.^(k-1)/(k-1)!

% put in correct taus 
tau_2 = (c.^2-D*(-ll).^2)/factorial(2) - A*c;
tau_3 = (c.^3-D*(-ll).^3)/factorial(3) - A*(c.^2)/factorial(2);
tau_4 = (c.^4-D*(-ll).^4)/factorial(4) - A*(c.^3)/factorial(3);
tau_5 = (c.^5-D*(-ll).^5)/factorial(5) - A*(c.^4)/factorial(4);
tau_6 = (c.^6-D*(-ll).^6)/factorial(6) - A*(c.^5)/factorial(5);
tau_7 = (c.^7-D*(-ll).^7)/factorial(7) - A*(c.^6)/factorial(6);
tau_8 = (c.^8-D*(-ll).^8)/factorial(8) - A*(c.^7)/factorial(7);

C=diag(c);

coneq(1)=sum(b)-(1-(-ll')*theta');

  

if p==2
    coneq(2)=b'*c-(1-(-ll').^2*theta')./2;
elseif p==3
  coneq(2)=b'*c-(1-(-ll').^2*theta')/2;
  coneq(3)=b'*c.^2-(1-(-ll').^(3)*theta')/3;
  coneq(4)=b'*tau_2;
elseif p==4
    coneq(2)=b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)=b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)=b'*tau_2;
    coneq(5)=b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(6)=b'*C*tau_2;
    coneq(7)=b'*A*tau_2;
    coneq(8)=b'*tau_3;
elseif p==5          %p>=5 tau_2=0
    coneq(2)=b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)=b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)=b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(5)=b'*tau_3;
    coneq(6)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(7)=b'*A*tau_3;
    coneq(8)=b'*tau_4;
    coneq(9)=b'*C*tau_3;
    coneq=[coneq tau_2'];
elseif p==6 
    coneq(2)=b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)=b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)=b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(5)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(6)=b'*tau_3;
    coneq(7)=b'*A*tau_3;
    coneq(8)=b'*tau_4;
    coneq(9)=b'*C*tau_3;
    coneq(10)=b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(11)=b'*A^2*tau_3;
    coneq(12)=b'*A*tau_4;
    coneq(13)=b'*A*C*tau_3;
    coneq(14)=b'*tau_5;
    coneq(15)=b'*C*A*tau_3;
    coneq(16)=b'*C*tau_4;
    coneq(17)=b'*C^2*tau_3;
    coneq=[coneq tau_2'];
elseif p==7   %p>=7 tau_2=0 & tau_3=0
    coneq(2)=b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)=b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)=b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(5)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(6)=b'*tau_4;
    coneq(7)=b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(8)=b'*A*tau_4;
    coneq(9)=b'*tau_5;
    coneq(10)=b'*C*tau_4;
    coneq(11)=b'*c.^6-(1-(-ll').^7*theta')/7;
    coneq(12)=b'*A^2*tau_4;
    coneq(13)=b'*A*tau_5;
    coneq(14)=b'*A*C*tau_4;
    coneq(15)=b'*tau_6;
    coneq(16)=b'*C*A*tau_4;
    coneq(17)=b'*C*tau_5;
    coneq(18)=b'*C^2*tau_4;
    coneq=[coneq tau_2' tau_3'];
elseif p==8
    coneq(2)=b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)=b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)=b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(5)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(6)=b'*tau_4;
    coneq(7)=b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(8)=b'*A*tau_4;
    coneq(9)=b'*tau_5;
    coneq(10)=b'*C*tau_4;
    coneq(11)=b'*c.^6-(1-(-ll').^7*theta')/7;
    coneq(12)=b'*A^2*tau_4;
    coneq(13)=b'*A*tau_5;
    coneq(14)=b'*A*C*tau_4;
    coneq(15)=b'*tau_6;
    coneq(16)=b'*C*A*tau_4;
    coneq(17)=b'*C*tau_5;
    coneq(18)=b'*C^2*tau_4;
    coneq(19) = b'*c.^7-(1-(-ll').^8*theta')/8;
    coneq(21) = b'*A^3*tau_4;
    coneq(22) = b'*A^2*tau_5;
    coneq(23) = b'*A^2*C*tau_4;
    coneq(24) = b'*A*tau_6;
    coneq(25) = b'*A*C*A*tau_4;
    coneq(26) = b'*A*C*tau_5;
    coneq(27) = b'*A*C^2*tau_4;
    coneq(28) = b'*tau_7;
    coneq(29) = b'*C*A^2*tau_4;
    coneq(30) = b'*C*A*tau_5;
    coneq(31) = b'*C*A*C*tau_4;
    coneq(32) = b'*C*tau_6;
    coneq(33) = b'*C^2*A*tau_4;
    coneq(34) = b'*C^2*tau_5;
    coneq(35) = b'*C^3*tau_4;

    coneq=[coneq tau_2' tau_3'];
elseif p>8 %p>=7 tau_2=0  tau_3=0 tau_4=0
    disp('Order conditions for p>8 are not coded up yet');
end
