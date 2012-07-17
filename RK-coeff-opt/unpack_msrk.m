function [A,Ahat,b,bhat,D,theta] =  unpack_msrk(X,s,k,class)
% function [A,Ahat,b,bhat,D,theta] =  unpack_msrk(X,s,k,class)
%
% Extract the coefficient arrays from the optimization vector

switch class 
    % TODO: clean up the following code
    case {'emsrk1', 'emsrk2'}
        count1= (s^2-s)/2;
        A=zeros(s-1,s);
        count=1;
        for i = 1:s-1
            A(i,1:i)=X(count:count+i-1);
            count= count+i;
        end
        A=[zeros(1,s);A];
    case {'imsrk1', 'imsrk2'}
        count1=s^2;
        A=reshape(X(1:s^2),s,s);
    case {'dimsrk1', 'dimsrk2'}
        count1=(s^2+s)/2;
        count=0;
        A=zeros(s,s);
        for i=1:s
            for j=1:i
                count=count+1;
                A(i,j)=X(count);
            end
        end
end

switch class
    case {'emsrk2', 'imsrk2', 'dimsrk2'}  % Type two methods
        count2= count1 + (s-1)*(k-1);  %keeps track the size of each vector
        count3= count2+s;              % Count keeps track of index while the 
        count4= count3 + (k-1);        %the addition is the number of unknowns in that vector 
        count5 = count4 +(s-1)*(k-1);
        count6= count5+k-1;
    
        Ahat=zeros(1,k-1); 
        Ahat= [Ahat; reshape(X(count1+1:count2),s-1,k-1)];
        b= X(count2+1:count3);
        bhat= X(count3+1:count4);
        D=reshape(X(count4+1:count5),s-1,k-1);
        D(:,end+1) = 1- sum(D,2);
        theta = X(count5+1:count6);
        theta(end+1)= 1-sum(theta);

    case {'emsrk1', 'imsrk1', 'dimsrk1'}  % Type one methods
        Ahat= [];
        bhat=[];
    
        D=reshape(X(count1+1:(count1+s*k-s)),s,k-1);
        b=reshape(X(count1+s*k+1-s:count1+s*k),1,s);
        theta=reshape(X((count1+s*k+1):end-1),1,k-1);
		  
        D(:,k)=1-sum(D,2);
        theta(k)=1-sum(theta);
end
