function coneq=oc_albrecht(A,b,c,p)
% function coneq=oc_albrecht(A,b,c,p)
%
% Order conditions for SSP RK Methods.
%
% This version is based on Albrecht's approach

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

%min_stage_order = floor((p-1)/2.);
min_stage_order = 1;

C=diag(c);

%=====================================================
% Stage order conditions

tau_2 = (c.^2)/factorial(2) - A*c;
tau_3 = (c.^3)/factorial(3) - A*(c.^2)/factorial(2);
tau_4 = (c.^4)/factorial(4) - A*(c.^3)/factorial(3);
tau_5 = (c.^5)/factorial(5) - A*(c.^4)/factorial(4);
tau_6 = (c.^6)/factorial(6) - A*(c.^5)/factorial(5);
tau_7 = (c.^7)/factorial(7) - A*(c.^6)/factorial(6);

%=====================================================
% Order conditions

% For most classes, sum(b)=1 is already imposed as a linear constraint
% But for low-storage methods it is not, and it does no harm
% to add it here in both cases.
coneq(1)=sum(b)-1.;

if min_stage_order<2

    if p>=2
        coneq(2)=b'*c-1./2;
    end

    if p>=3
      coneq(3)=b'*c.^2-1./3;
      coneq(4)=b'*tau_2;
    end

    if p>=4 
      coneq(5)=b'*c.^3-1./4;
      coneq(6)=b'*C*tau_2;
      coneq(7)=b'*A*tau_2;
      coneq(8)=b'*tau_3;
    end

    if p>=5 
      coneq(9)=b'*c.^4-1./5;
      coneq(10)=b'*A*C*tau_2;
      coneq(11)=b'*A^2*tau_2;
      coneq(12)=b'*A*tau_3;
      coneq(13)=b'*tau_4;
      coneq(14)=b'*C*A*tau_2;
      coneq(15)=b'*C*tau_3;
      coneq(16)=b'*C^2*tau_2;
      coneq(17)=b'*tau_2.^2;
    end

    if p>=6 
        coneq(18)=b'*c.^5-1./6;
        coneq(19)=b'*A^2*C*tau_2;
        coneq(20)=b'*A^3*tau_2;
        coneq(21)=b'*A^2*tau_3;
        coneq(22)=b'*A*tau_4;
        coneq(23)=b'*A*C*A*tau_2;
        coneq(24)=b'*A*C*tau_3;
        coneq(25)=b'*A*C^2*tau_2;
        coneq(26)=b'*A*(tau_2.^2);
        coneq(27)=b'*tau_5;
        coneq(28)=b'*C*A*C*tau_2;
        coneq(29)=b'*C*A^2*tau_2;
        coneq(30)=b'*C*A*tau_3;
        coneq(31)=b'*C*tau_4;
        coneq(32)=b'*C^2*A*tau_2;
        coneq(33)=b'*C^2*tau_3;
        coneq(34)=b'*C^3*tau_2;
        coneq(35)=b'*C*(tau_2.^2);
        coneq(36)=b'*(tau_2.*tau_3);
        coneq(37)=b'*(tau_2.*(A*tau_2));
    end

    if p>6
        disp('Albrecht-form order conditions for p>6 not coded yet')
    end

elseif min_stage_order<3 % min stage order == 2
      coneq(2)=b'*c-1./2;
      coneq(3)=b'*c.^2-1./3;
      coneq(4)=b'*c.^3-1./4;
      coneq(5)=b'*tau_3;
      coneq(6)=b'*c.^4-1./5;
      coneq(7)=b'*A*tau_3;
      coneq(8)=b'*tau_4;
      coneq(9)=b'*C*tau_3;

    if p>=6 
        coneq(10)=b'*c.^5-1./6;
        coneq(11)=b'*A^2*tau_3;
        coneq(12)=b'*A*tau_4;
        coneq(13)=b'*A*C*tau_3;
        coneq(14)=b'*tau_5;
        coneq(15)=b'*C*A*tau_3;
        coneq(16)=b'*C*tau_4;
        coneq(17)=b'*C^2*tau_3;
    end
    coneq=[coneq tau_2'];

    if p>6
        disp('Albrecht-form order conditions for p>6 not coded yet')
    end

elseif min_stage_order==3
      coneq(2)=b'*c   -1./2;
      coneq(3)=b'*c.^2-1./3;
      coneq(4)=b'*c.^3-1./4;
      coneq(5)=b'*c.^4-1./5;
      coneq(6)=b'*tau_4;
      coneq(7)=b'*c.^5-1./6;
      coneq(8)=b'*A*tau_4;
      coneq(9)=b'*tau_5;
      coneq(10)=b'*C*tau_4;
      coneq(11)=b'*c.^6-1./7;
      coneq(12)=b'*A^2*tau_4;
      coneq(13)=b'*A*tau_5;
      coneq(14)=b'*A*C*tau_4;
      coneq(15)=b'*tau_6;
      coneq(16)=b'*C*A*tau_4;
      coneq(17)=b'*C*tau_5;
      coneq(18)=b'*C^2*tau_4;

    if p>=8
        coneq(19) = b'*c.^7-1./8;
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
    end

    if p==9
        coneq(36) = b'*tau_8;
        coneq(37) = b'*A*tau_7;
        coneq(38) = b'*A^2*tau_6;
        coneq(39) = b'*A^3*tau_5;
        coneq(40) = b'*A^4*tau_4;
        coneq(41) = b'*A^3*C*tau_4;
        coneq(42) = b'*A^2*C*tau_5;
        coneq(43) = b'*A^2*C*A*tau_4;
        coneq(44) = b'*A^2*C^2*tau_4;
        coneq(45) = b'*A*C*tau_6;
        coneq(46) = b'*A*C*A*tau_5;
        coneq(47) = b'*A*C*A^2*tau_4;
        coneq(48) = b'*A*C*A*C*tau_4;
        coneq(49) = b'*A*C^2*tau_5;
        coneq(50) = b'*A*C^2*A*tau_4;
        coneq(51) = b'*A*C^3*tau_4;
        coneq(52) = b'*C*tau_7;
        coneq(53) = b'*C*A*tau_6;
        coneq(54) = b'*C*A*A*tau_5;
        coneq(55) = b'*C*A*A*A*tau_4;
        coneq(56) = b'*C*A*A*C*tau_4;
        coneq(57) = b'*C*A*C*tau_5;
        coneq(58) = b'*C*A*C*A*tau_4;
        coneq(59) = b'*C*A*C^2*tau_4;
        coneq(60) = b'*C^2*tau_6;
        coneq(61) = b'*C^2*A*tau_5;
        coneq(62) = b'*C^2*A^2*tau_4;
        coneq(63) = b'*C^2*A*C*tau_4;
        coneq(64) = b'*C^3*tau_5;
        coneq(65) = b'*C^3*A*tau_4;
        coneq(66) = b'*C^4*tau_4;
        coneq(67) = b'*tau_4.^2;
    end

    coneq=[coneq tau_2' tau_3'];
    if p>9
        disp('Albrecht-form order conditions for p>6 not coded yet')
    end
else disp('Order conditions for stage_order>3 are not coded up yet');
end
end
