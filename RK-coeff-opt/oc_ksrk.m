function coneq= oc_ksrk(A,b,D,theta,p)
% function coneq= oc_ksrk(A,b,D,theta,p)
% Order conditions for multistep-RK methods.
%
% ..warning:: 
%
%         Here we assume a certain minimum stage order,
%         which is necessarily true for methods with
%         strictly positive abscissae (b>0).
%         This assumption dramatically reduces the
%         number of order conditions that must be
%         considered for high-order methods.
%         For methods that do not satisfy b>0, this
%         assumption may be unnecessarily restrictive.


k=length(theta);
ll=[k-1:-1:0]';
c=A*ones(length(b),1)-D*ll;
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

% First order condition; this should really be
% implemented with the linear constraints
coneq(1)=sum(b)-(1-(-ll')*theta');

if p>=2
    coneq(2)=b'*c   -(1-(-ll').^2*theta')/2;
if p>=3
    coneq(3)=b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)=b'*tau_2;
if p>=4
    coneq(5)=b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(6)=b'*C*tau_2;
    coneq(7)=b'*A*tau_2;
    coneq(8)=b'*tau_3;
if p==5          % for p>=5 we assume tau_2=0
    coneq(9)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(10)=b'*A*tau_3;
    coneq(11)=b'*tau_4;
    coneq(12)=b'*C*tau_3;
    coneq=[coneq tau_2'];
elseif p==6 
    coneq(9)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(10)=b'*A*tau_3;
    coneq(11)=b'*tau_4;
    coneq(12)=b'*C*tau_3;
    coneq(13)=b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(14)=b'*A^2*tau_3;
    coneq(15)=b'*A*tau_4;
    coneq(16)=b'*A*C*tau_3;
    coneq(17)=b'*tau_5;
    coneq(18)=b'*C*A*tau_3;
    coneq(19)=b'*C*tau_4;
    coneq(20)=b'*C^2*tau_3;
    coneq=[coneq tau_2'];
elseif p==7   %for p>=7 we assume tau_2=0 & tau_3=0
    coneq(9)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(10)=b'*tau_4;
    coneq(11)=b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(12)=b'*A*tau_4;
    coneq(13)=b'*tau_5;
    coneq(14)=b'*C*tau_4;
    coneq(15)=b'*c.^6-(1-(-ll').^7*theta')/7;
    coneq(16)=b'*A^2*tau_4;
    coneq(17)=b'*A*tau_5;
    coneq(18)=b'*A*C*tau_4;
    coneq(19)=b'*tau_6;
    coneq(20)=b'*C*A*tau_4;
    coneq(21)=b'*C*tau_5;
    coneq(22)=b'*C^2*tau_4;
    coneq=[coneq tau_2' tau_3'];
elseif p==8
    coneq(9)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(10)=b'*tau_4;
    coneq(11)=b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(12)=b'*A*tau_4;
    coneq(13)=b'*tau_5;
    coneq(14)=b'*C*tau_4;
    coneq(15)=b'*c.^6-(1-(-ll').^7*theta')/7;
    coneq(16)=b'*A^2*tau_4;
    coneq(17)=b'*A*tau_5;
    coneq(18)=b'*A*C*tau_4;
    coneq(19)=b'*tau_6;
    coneq(20)=b'*C*A*tau_4;
    coneq(21)=b'*C*tau_5;
    coneq(22)=b'*C^2*tau_4;
    coneq(23) = b'*c.^7-(1-(-ll').^8*theta')/8;
    coneq(24) = b'*A^3*tau_4;
    coneq(25) = b'*A^2*tau_5;
    coneq(26) = b'*A^2*C*tau_4;
    coneq(27) = b'*A*tau_6;
    coneq(28) = b'*A*C*A*tau_4;
    coneq(29) = b'*A*C*tau_5;
    coneq(30) = b'*A*C^2*tau_4;
    coneq(31) = b'*tau_7;
    coneq(32) = b'*C*A^2*tau_4;
    coneq(33) = b'*C*A*tau_5;
    coneq(34) = b'*C*A*C*tau_4;
    coneq(35) = b'*C*tau_6;
    coneq(36) = b'*C^2*A*tau_4;
    coneq(37) = b'*C^2*tau_5;
    coneq(38) = b'*C^3*tau_4;

    coneq=[coneq tau_2' tau_3'];
elseif p>8 %for p>=9 we assume tau_2=0  tau_3=0 tau_4=0
    disp('Order conditions for p>8 are not coded up yet');
end
