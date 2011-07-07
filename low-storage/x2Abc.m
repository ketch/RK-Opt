function [A,b,bhat,c,alpha,beta,gamma1,gamma2,gamma3,delta]=x2Abc(regs,class,x)
%function [A,b,bhat,c,alpha,beta]=x2Abc(regs,class,x)
%Convert vector x of low-storage coefficients to Butcher array

gamma3=[];

switch regs
  case 2 %Two-register methods

    switch class
      case 'plain'
        % n = 3s - 3 free parameters
        s=length(x)/3+1; % # of stages
        %x=[gamma_{32} ... gamma_{s+1,2} 
        %            delta2 ... delta_{s-1} beta21 ... beta_{s+1,s}]
        delta =[1 x(s:2*s-3)' 0];
      case 'star'
        s=(length(x))/3+1;
        %x=[gamma_{32} ... gamma_{s+1,2} 
        %       delta2 ... delta_{s-1} beta21 ... beta_{s+1,s}]
        delta =[1 zeros(size(x(s:2*s-3)')) 0];
      case 'emb'
        % n = 3s - 1 free parameters
        s=(length(x)+1)/3; % # of stages
        %x=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1} 
                    % beta21 ... beta_{s+1,s}   delta_s delta_{s+1}]
        delta =[1 x(s:2*s-3)' x(end-1) x(end)];
    end

    gamma2=[0 1 x(1:s-1)'];

    for i=2:s+1 gamma1(i)=1-gamma2(i)*sum(delta(1:i-1)); end

    alpha=zeros(s+1,s); beta=zeros(s+1,s);
    for i=1:s beta(i+1,i)=x(2*s-3+i); end

    for i=2:s
      beta( i+1,i-1) = -beta(i,i-1)*gamma2(i+1)/gamma2(i);
      alpha(i+1,i-1) = -gamma1(i)  *gamma2(i+1)/gamma2(i);
      alpha(i+1,i  ) = 1 - alpha(i+1,i-1);
    end

  case 3 %Three-register methods

    switch class
      case 'star'
        % n = 4s - 6 free parameters
        s=(length(x)+6)/4; % # of stages
        %for 3S* methods witohout an embedded method
        %x=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1} 
        %    beta21 ... beta_{s+1,s} gamma_{53} ... gamma_{s+1,3} ]
        delta =[1 x(s:2*s-3)' 0];
      case 'staremb'
        % n = 4s - 3 free parameters
        s=(length(x)+3)/4; % # of stages
        %for 3S* methods with an embedded method
        %x=[gamma_{32} ... gamma_{s+1,2} delta2 ... delta_{s-1} 
        %      beta21 ... beta_{s+1,s}
        %     gamma_{53} ... gamma_{s+1,3} delta_s delta_{s+1} delta_{s+2} ]
        delta =[1 x(s:2*s-3)' x(end-2:end)'];
    end

    alpha=zeros(s+1,s); beta=zeros(s+1,s);
    for i=1:s beta(i+1,i)=x(2*s-3+i); end

    gamma2=[0 1 x(1:s-1)'];
    gamma3=[0 0 0 0 x(3*s-2:4*s-6)'];

    for i=2:s+1 gamma1(i)=1-gamma3(i)-gamma2(i)*sum(delta(1:i-1)); end
    %alpha(2,1)=gamma1(2)+gamma2(2)+gamma3(2); %Should be 1

    for i=2:s
      beta( i+1,i-1) = -beta(i,i-1)*gamma2(i+1)/gamma2(i);
      alpha(i+1,1  ) = gamma3(i+1)-gamma3(i)*gamma2(i+1)/gamma2(i);
      alpha(i+1,i-1) = -gamma1(i)  *gamma2(i+1)/gamma2(i);
      alpha(i+1,i  ) = 1 - alpha(i+1,i-1) - alpha(i+1,1);
    end

end

%Convert Shu-Osher to Butcher
[A,b,c]=shuosher2butcher(alpha,beta);

%Now do embedded method
if strcmp(class,'emb') && regs==2
    bhat = (delta*[A;b']/sum(delta))';
elseif strcmp(class,'staremb') && regs==3
  bhat = (delta(1:s+1)*[A;b']/sum(delta))';
else
  bhat = [];
end
