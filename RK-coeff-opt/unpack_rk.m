function [A,b,c]=unpack_rk(X,s,class)
% function [A,b,c]=unpack_rk(X,s,class)
%
% Extracts the coefficient arrays from the optimization vector.
%
% The coefficients are tored in a single vector x as::
%
%       x=[A b' c']
%
% A is stored row-by-row.

A=zeros(s);
%Generate Butcher coefficients (for RK) or Shu-Osher coefficients (for lsRK)
switch class 
    case 'erk'
        c=[0 X(1:s-1)]'; b=X(s:2*s-1)'; 
        for i=1:s
            A(i,1:i-1)=X(2*s+(i-2)*(i-1)/2:2*s-1+i*(i-1)/2);
        end
    case 'irk'
        c=X(1:s)'; b=X(s+1:2*s)';
        A=reshape(X(2*s+1:2*s+s^2),s,s)';
    case 'irk5'
        c=[0 X(1:s-1)]'; b=X(s:2*s-1)';
        for i=2:s
            A(i,:)=X(i*s:(i+1)*s-1);
        end
    case 'dirk'
        c=X(1:s)'; b=X(s+1:2*s)';
        for i=1:s
            A(i,1:i)=X(2*s+1+i*(i-1)/2:2*s+i*(i+1)/2);
        end
    case 'dirk5'
        c=[0 X(1:s-1)]'; b=X(s:2*s-1)'; 
        for i=2:s
            A(i,1:i)=X(2*s-1+i*(i-1)/2:2*s-2+i*(i+1)/2);
        end
    case 'sdirk'
        c=X(1:s)'; b=X(s+1:2*s)';
        A(1,1)=X(2*s+1);
        for i=2:s
            A(i,1:i-1)=X(2*s+2+(i-2)*(i-1)/2:2*s+1+i*(i-1)/2);
            A(i,i)=X(2*s+1);
        end
    case '2S'
        % n = 3s - 3 free parameters
        s=length(X)/3+1; % # of stages
        %X=[gamma_{32} ... gamma_{s+1,2} 
        %            delta2 ... delta_{s-1} beta21 ... beta_{s+1,s}]
        delta =[1 X(s:2*s-3) 0];
    case '2Sstar'
        s=(length(X))/3+1;
        %X=[gamma_{32} ... gamma_{s+1,2} 
        %       delta2 ... delta_{s-1} beta21 ... beta_{s+1,s}]
        delta =[1 zeros(size(X(s:2*s-3))) 0];
    case '2Semb'
        % n = 3s - 1 free parameters
        s=(length(X)+1)/3; % # of stages
        %X=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1} 
                    % beta21 ... beta_{s+1,s}   delta_s delta_{s+1}]
        delta =[1 X(s:2*s-3) X(end-1) X(end)];
    case '3Sstar'
        % n = 4s - 6 free parameters
        s=(length(X)+6)/4; % # of stages
        assert(s>=3, '3S* methods require at least 3 stages');
        %for 3S* methods witohout an embedded method
        %X=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1} 
        %    beta21 ... beta_{s+1,s} gamma_{53} ... gamma_{s+1,3} ]
        delta =[1 X(s:2*s-3) 0];
    case '3Sstaremb'
        % n = 4s - 3 free parameters
        s=(length(X)+3)/4; % # of stages
        assert(s>=3, '3S* methods require at least 3 stages');
        %for 3S* methods with an embedded method
        %X=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1} 
        %      beta21 ... beta_{s+1,s}
        %     gamma_{53} ... gamma_{s+1,3} delta_s delta_{s+1} delta_{s+2} ]
        delta =[1 X(s:2*s-3) X(end-2:end)];
end


% Now generate Butcher coefficients for low-storage methods
if strcmp(class(1:2),'2S')
    gamma2=[0 1 X(1:s-1)];

    for i=2:s+1
        gamma1(i)=1-gamma2(i)*sum(delta(1:i-1)); 
    end

    alpha=zeros(s+1,s); beta=zeros(s+1,s);
    for i=1:s 
        beta(i+1,i)=X(2*s-3+i); 
    end

    for i=2:s
      beta( i+1,i-1) = -beta(i,i-1)*gamma2(i+1)/gamma2(i);
      alpha(i+1,i-1) = -gamma1(i)  *gamma2(i+1)/gamma2(i);
      alpha(i+1,i  ) = 1 - alpha(i+1,i-1);
    end

    %Convert Shu-Osher to Butcher
    [A,b,c]=shuosher2butcher(alpha,beta);

elseif strcmp(class(1:2),'3S')

    alpha=zeros(s+1,s); beta=zeros(s+1,s);
    for i=1:s 
        beta(i+1,i)=X(2*s-3+i); 
    end

    gamma2=[0 1 X(1:s-1)];
    gamma3=[0 0 0 0 X(3*s-2:4*s-6)];

    for i=2:s+1 
        gamma1(i)=1-gamma3(i)-gamma2(i)*sum(delta(1:i-1)); 
    end

    for i=2:s
      beta( i+1,i-1) = -beta(i,i-1)*gamma2(i+1)/gamma2(i);
      alpha(i+1,1  ) = gamma3(i+1)-gamma3(i)*gamma2(i+1)/gamma2(i);
      alpha(i+1,i-1) = -gamma1(i)  *gamma2(i+1)/gamma2(i);
      alpha(i+1,i  ) = 1 - alpha(i+1,i-1) - alpha(i+1,1);
    end

    %Convert Shu-Osher to Butcher
    [A,b,c]=shuosher2butcher(alpha,beta);
end
