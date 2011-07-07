function [A,b,c]=unpack_rk(X,s,class)
%function [A,b,c]=unpack_rk(X,s,class)
%
% Extracts the coefficient arrays from the optimization vector
% Stored in a single vector x as:
% x=[A b' c']
% A is stored row-by-row

A=zeros(s);
switch class 
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
  case 'erk'
    c=[0 X(1:s-1)]'; b=X(s:2*s-1)'; 
    for i=1:s
      A(i,1:i-1)=X(2*s+(i-2)*(i-1)/2:2*s-1+i*(i-1)/2);
    end
end
