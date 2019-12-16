function [A,Ahat,b,bhat,c,chat,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,class)
% function [A,b,bhat,c,alpha,beta,gamma1,gamma2,gamma3,delta]=unpack_lsrk(X,class)
%
% Extracts the coefficient arrays from the optimization vector.
%
% This function also returns the low-storage coefficients.

gamma3=[];
switch class
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
    %for 3S* methods witohout an embedded method
    %X=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1}
    %    beta21 ... beta_{s+1,s} gamma_{53} ... gamma_{s+1,3} ]
    delta =[1 X(s:2*s-3) 0];
  case '3Sstaremb'
    % n = 4s - 3 free parameters
    s=(length(X)+3)/4; % # of stages
    %for 3S* methods with an embedded method
    %X=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1}
    %      beta21 ... beta_{s+1,s}
    %     gamma_{53} ... gamma_{s+1,3} delta_s delta_{s+1} delta_{s+2} ]
    delta =[1 X(s:2*s-3) X(end-2:end)];
  case '3SstarembFSAL'
    % n = 4s - 2 free parameters
    s=(length(X)+2)/4; % # of stages
    assert(s>=3, '3S* methods require at least 3 stages');
    %for 3S* methods with an embedded method
    %X=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1}
    %      beta21 ... beta_{s+1,s}
    %     gamma_{53} ... gamma_{s+1,3} delta_s delta_{s+1} delta_{s+2} ]
    delta =[1 X(s:2*s-3) X(end-3:end-1)];
    betahat = X(end);
end

if strcmp(class(1:2),'2S')
    gamma2=[0 1 X(1:s-1)];

    gamma1 = zeros(1,s+1);
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

    gamma1 = zeros(1,s+1);
    for i=2:s+1
        gamma1(i)=1-gamma3(i)-gamma2(i)*sum(delta(1:i-1));
    end
    %alpha(2,1)=gamma1(2)+gamma2(2)+gamma3(2); %Should be 1

    for i=2:s
        beta( i+1,i-1) = -beta(i,i-1)*gamma2(i+1)/gamma2(i);
        alpha(i+1,1  ) = gamma3(i+1)-gamma3(i)*gamma2(i+1)/gamma2(i);
        alpha(i+1,i-1) = -gamma1(i)  *gamma2(i+1)/gamma2(i);
        alpha(i+1,i  ) = 1 - alpha(i+1,i-1) - alpha(i+1,1);
    end

    %Convert Shu-Osher to Butcher
    [A,b,c]=shuosher2butcher(alpha,beta);
end

%Now do embedded method
Ahat = A;
chat = c;
if strcmp(class,'2Semb')
    bhat = (delta*[A;b']/sum(delta))';
elseif strcmp(class,'3Sstaremb')
  bhat = (delta(1:s+1)*[A;b']/sum(delta))';
elseif strcmp(class,'3SstarembFSAL')
  bhat = (delta(1:s+1)*[A;b']/sum(delta))';
  bhat(end+1) = betahat;
  Ahat = A;
  Ahat(end+1,:) = b';
  Ahat(:,end+1) = 0;
  chat(end+1) = 1;
else
  bhat = [];
end
