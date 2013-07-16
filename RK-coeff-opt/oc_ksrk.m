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

tau_2 = (c.^2-D*(-ll).^2)/factorial(2) - A*c;
tau_3 = (c.^3-D*(-ll).^3)/factorial(3) - A*(c.^2)/factorial(2);
tau_4 = (c.^4-D*(-ll).^4)/factorial(4) - A*(c.^3)/factorial(3);
tau_5 = (c.^5-D*(-ll).^5)/factorial(5) - A*(c.^4)/factorial(4);
tau_6 = (c.^6-D*(-ll).^6)/factorial(6) - A*(c.^5)/factorial(5);
tau_7 = (c.^7-D*(-ll).^7)/factorial(7) - A*(c.^6)/factorial(6);
tau_8 = (c.^8-D*(-ll).^8)/factorial(8) - A*(c.^7)/factorial(7);
tau_9 = (c.^9-D*(-ll).^9)/factorial(9) - A*(c.^8)/factorial(8);
tau_10 = (c.^10-D*(-ll).^10)/factorial(10) - A*(c.^9)/factorial(9);
tau_11 = (c.^11-D*(-ll).^11)/factorial(11) - A*(c.^10)/factorial(10);

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

elseif p==8 % Assume stage order 3
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

elseif p==9 % Assume stage order 4
    coneq(2)=b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)=b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)=b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(5)=b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(6)=b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(7)=b'*tau_5;
    coneq(8)=b'*c.^6-(1-(-ll').^7*theta')/7;
    coneq(9)=b'*A*tau_5;
    coneq(10)=b'*tau_6;
    coneq(11)=b'*C*tau_5;
    coneq(12) = b'*c.^7-(1-(-ll').^8*theta')/8;
    coneq(13) = b'*A^2*tau_5;
    coneq(14) = b'*A*tau_6;
    coneq(15) = b'*A*C*tau_5;
    coneq(16) = b'*tau_7;
    coneq(17) = b'*C*A*tau_5;
    coneq(18) = b'*C*tau_6;
    coneq(19) = b'*C^2*tau_5;
    coneq(20) = b'*c.^8-(1-(-ll').^9*theta')/9;
    coneq(21) = b'*C*A^2*tau_5;
    coneq(22) = b'*C*A*tau_6;
    coneq(23) = b'*C*A*C*tau_5;
    coneq(24) = b'*C*tau_7;
    coneq(25) = b'*C^2*A*tau_5;
    coneq(26) = b'*C^2*tau_6;
    coneq(27) = b'*C^3*tau_5;
    coneq(28) = b'*A^2*tau_5;
    coneq(29) = b'*A^3*tau_5;
    coneq(30) = b'*A^2*tau_6;
    coneq(31) = b'*A^2*C*tau_5;
    coneq(32) = b'*A*tau_7;
    coneq(33) = b'*A*C*A*tau_5;
    coneq(34) = b'*A*C*tau_6;
    coneq(35) = b'*A*C^2*tau_5;
    coneq(36) = b'*tau_8;

    coneq=[coneq tau_2' tau_3' tau_4'];

elseif p==10  % Assume stage order 4
    coneq(2)  = b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)  = b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)  = b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(5)  = b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(6)  = b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(7)  = b'*c.^6-(1-(-ll').^7*theta')/7;
    coneq(8)  = b'*c.^7-(1-(-ll').^8*theta')/8;	
    coneq(9)  = b'*c.^8-(1-(-ll').^9*theta')/9;    
    coneq(10) = b'*A*tau_5;
    coneq(11) = b'*tau_6;
    coneq(12) = b'*C*tau_5;
    coneq(13) = b'*tau_5;
    coneq(14) = b'*A^2*tau_5;
    coneq(15) = b'*A*tau_6;
    coneq(16) = b'*A*C*tau_5;
    coneq(17) = b'*tau_7;
    coneq(18) = b'*C*A*tau_5;
    coneq(19) = b'*C*tau_6;
    coneq(20) = b'*C^2*tau_5;
    coneq(21) = b'*C*A^2*tau_5;
    coneq(22) = b'*C*A*tau_6;
    coneq(23) = b'*C*A*C*tau_5;
    coneq(24) = b'*C*tau_7;
    coneq(25) = b'*C^2*A*tau_5;
    coneq(26) = b'*C^2*tau_6;
    coneq(27) = b'*C^3*tau_5;
    coneq(28) = b'*A^2*tau_5;
    coneq(29) = b'*A^3*tau_5;
    coneq(30) = b'*A^2*tau_6;
    coneq(31) = b'*A^2*C*tau_5;
    coneq(32) = b'*A*tau_7;
    coneq(33) = b'*A*C*A*tau_5;
    coneq(34) = b'*A*C*tau_6;
    coneq(35) = b'*A*C^2*tau_5;
    coneq(36) = b'*tau_8;
    coneq(37) = b'*c.^9-(1-(-ll').^10*theta')/10; 
    coneq(38) = b'*A^2*tau_5;
    coneq(39) = b'*A*tau_6;
    coneq(40) = b'*A*C*tau_5;
    coneq(41) = b'*A*tau_5;
    coneq(42) = b'*A^3*tau_5;
    coneq(43) = b'*A^2*tau_6;
    coneq(44) = b'*A^2*C*tau_5;
    coneq(45) = b'*A*tau_7;
    coneq(46) = b'*A*C*A*tau_5;
    coneq(47) = b'*A*C*tau_6;
    coneq(48) = b'*A*C^2*tau_5;
    coneq(49) = b'*A*C*A^2*tau_5;
    coneq(50) = b'*A*C*A*tau_6;
    coneq(51) = b'*A*C*A*C*tau_5;
    coneq(52) = b'*A*C*tau_7;
    coneq(53) = b'*A*C^2*A*tau_5;
    coneq(54) = b'*A*C^2*tau_6;
    coneq(55) = b'*A*C^3*tau_5;
    coneq(56) = b'*A^3*tau_5;
    coneq(57) = b'*A^4*tau_5;
    coneq(58) = b'*A^3*tau_6;
    coneq(59) = b'*A^3*C*tau_5;
    coneq(60) = b'*A^2*tau_7;
    coneq(61) = b'*A^2*C*A*tau_5;
    coneq(62) = b'*A^2*C*tau_6;
    coneq(63) = b'*A^3*tau_5;
    coneq(64) = b'*A*tau_8;
    coneq(65) = b'*C*A*tau_5;
    coneq(66) = b'*C*tau_6;
    coneq(67) = b'*C^2*tau_5;
    coneq(68) = b'*C*tau_5;
    coneq(69) = b'*C*A^2*tau_5;
    coneq(70) = b'*C*A*tau_6;
    coneq(71) = b'*C*A*C*tau_5;
    coneq(72) = b'*C*tau_7;
    coneq(73) = b'*C^2*A*tau_5;
    coneq(74) = b'*C^2*tau_6;
    coneq(75) = b'*C^3*tau_5;
    coneq(76) = b'*C^2*A^2*tau_5;
    coneq(77) = b'*C^2*A*tau_6;
    coneq(78) = b'*C^2*A*C*tau_5;
    coneq(79) = b'*C^2*tau_7;
    coneq(80) = b'*C^3*A*tau_5;
    coneq(81) = b'*C^3*tau_6;
    coneq(82) = b'*C^4*tau_5;
    coneq(83) = b'*C*A^2*tau_5;
    coneq(84) = b'*C*A^3*tau_5;
    coneq(85) = b'*C*A^2*tau_6;
    coneq(86) = b'*C*A^2*C*tau_5;
    coneq(87) = b'*C*A*tau_7;
    coneq(88) = b'*C*A*C*A*tau_5;
    coneq(89) = b'*C*A*C*tau_6;
    coneq(90) = b'*C*A*C^2*tau_5;
    coneq(91) = b'*C*tau_8;
    coneq(92) = b'*tau_9;

    coneq=[coneq tau_2' tau_3' tau_4'];
  
elseif p==11; % Assume stage order 5
    coneq(2)  = b'*c-(1-(-ll').^2*theta')/2;
    coneq(3)  = b'*c.^2-(1-(-ll').^3*theta')/3;
    coneq(4)  = b'*c.^3-(1-(-ll').^4*theta')/4;
    coneq(5)  = b'*c.^4-(1-(-ll').^5*theta')/5;
    coneq(6)  = b'*c.^5-(1-(-ll').^6*theta')/6;
    coneq(7)  = b'*c.^6-(1-(-ll').^7*theta')/7;
    coneq(8)  = b'*c.^7-(1-(-ll').^8*theta')/8;	
    coneq(9)  = b'*c.^8-(1-(-ll').^9*theta')/9;      
    coneq(10) = b'*c.^9-(1-(-ll').^10*theta')/10;
    coneq(11) = b'*c.^10-(1-(-ll').^11*theta')/11;

    coneq(12) = b'*tau_6;
    coneq(13) = b'*A*tau_6;
    coneq(14) = b'*tau_7;     
    coneq(15) = b'*C*tau_6;
    coneq(16) = b'*C*A*tau_6;
    coneq(17) = b'*C*tau_7;
    coneq(18) = b'*C^2*tau_6;
    coneq(19) = b'*A^2*tau_6;
    coneq(20) = b'*A*tau_7;       
    coneq(21) = b'*A*C*tau_6;       
    coneq(22) = b'*tau_8;
    coneq(23) = b'*A*tau_6;
    coneq(24) = b'*A^2*tau_6;
    coneq(25) = b'*A*tau_7;
    coneq(26) = b'*A*C*tau_6;
    coneq(27) = b'*A*C*A*tau_6;
    coneq(28) = b'*A*C*tau_7;
    coneq(29) = b'*A*C^2*tau_6;
    coneq(30) = b'*A^3*tau_6;
    coneq(31) = b'*A^2*tau_7;
    coneq(32) = b'*A^2*C*tau_6;
    coneq(33) = b'*A*tau_8;
    coneq(34) = b'*C*tau_6;
    coneq(35) = b'*C*A*tau_6;
    coneq(36) = b'*C*tau_7;
    coneq(37) = b'*C^2*tau_6;
    coneq(38) = b'*C^2*A*tau_6;
    coneq(39) = b'*C^2*tau_7;
    coneq(40) = b'*C^3*tau_6;
    coneq(41) = b'*C*A^2*tau_6;
    coneq(42) = b'*C*A*tau_7;
    coneq(43) = b'*C*A*C*tau_6;
    coneq(44) = b'*C*tau_8;
    coneq(45) = b'*tau_9;

    coneq(46) = b'*A*tau_6;
    coneq(47) = b'*A*A*tau_6;
    coneq(48) = b'*A*tau_7;     
    coneq(49) = b'*A*C*tau_6;
    coneq(50) = b'*A*C*A*tau_6;
    coneq(51) = b'*A*C*tau_7;
    coneq(52) = b'*A*C^2*tau_6;
    coneq(53) = b'*A*A^2*tau_6;
    coneq(54) = b'*A*A*tau_7;       
    coneq(55) = b'*A*A*C*tau_6;       
    coneq(56) = b'*A*tau_8;
    coneq(57) = b'*A*A*tau_6;
    coneq(58) = b'*A*A^2*tau_6;
    coneq(59) = b'*A*A*tau_7;
    coneq(60) = b'*A*A*C*tau_6;
    coneq(61) = b'*A*A*C*A*tau_6;
    coneq(62) = b'*A*A*C*tau_7;
    coneq(63) = b'*A*A*C^2*tau_6;
    coneq(64) = b'*A*A^3*tau_6;
    coneq(65) = b'*A*A^2*tau_7;
    coneq(66) = b'*A*A^2*C*tau_6;
    coneq(67) = b'*A*A*tau_8;
    coneq(68) = b'*A*C*tau_6;
    coneq(69) = b'*A*C*A*tau_6;
    coneq(70) = b'*A*C*tau_7;
    coneq(71) = b'*A*C^2*tau_6;
    coneq(72) = b'*A*C^2*A*tau_6;
    coneq(73) = b'*A*C^2*tau_7;
    coneq(74) = b'*A*C^3*tau_6;
    coneq(75) = b'*A*C*A^2*tau_6;
    coneq(76) = b'*A*C*A*tau_7;
    coneq(77) = b'*A*C*A*C*tau_6;
    coneq(78) = b'*A*C*tau_8;
    coneq(79) = b'*A*tau_9;     

    coneq(80) = b'*C*tau_6;
    coneq(81) = b'*C*A*tau_6;
    coneq(82) = b'*C*tau_7;     
    coneq(83) = b'*C*C*tau_6;
    coneq(84) = b'*C*C*A*tau_6;
    coneq(85) = b'*C*C*tau_7;
    coneq(86) = b'*C*C^2*tau_6;
    coneq(87) = b'*C*A^2*tau_6;
    coneq(88) = b'*C*A*tau_7;       
    coneq(89) = b'*C*A*C*tau_6;       
    coneq(90) = b'*C*tau_8;
    coneq(91) = b'*C*A*tau_6;
    coneq(92) = b'*C*A^2*tau_6;
    coneq(93) = b'*C*A*tau_7;
    coneq(94) = b'*C*A*C*tau_6;
    coneq(95) = b'*C*A*C*A*tau_6;
    coneq(96) = b'*C*A*C*tau_7;
    coneq(97) = b'*C*A*C^2*tau_6;
    coneq(98) = b'*C*A^3*tau_6;
    coneq(99) = b'*C*A^2*tau_7;
    coneq(100) = b'*C*A^2*C*tau_6;
    coneq(101) = b'*C*A*tau_8;
    coneq(102) = b'*C*C*tau_6;
    coneq(103) = b'*C*C*A*tau_6;
    coneq(104) = b'*C*C*tau_7;
    coneq(105) = b'*C*C^2*tau_6;
    coneq(106) = b'*C*C^2*A*tau_6;
    coneq(107) = b'*C*C^2*tau_7;
    coneq(108) = b'*C*C^3*tau_6;
    coneq(109) = b'*C*C*A^2*tau_6;
    coneq(110) = b'*C*C*A*tau_7;
    coneq(111) = b'*C*C*A*C*tau_6;
    coneq(112) = b'*C*C*tau_8;
    coneq(113) = b'*C*tau_9;
    coneq(114) = b'*tau_10;
        
    coneq=[coneq tau_2' tau_3' tau_4' tau_5'];

elseif p>11 
    disp('Order conditions for p>11 are not coded up yet');
end
