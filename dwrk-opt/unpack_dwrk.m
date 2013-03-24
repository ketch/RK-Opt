function [A,b,c,At,bt,ct]=unpack_dwrk(X,s,class)
% function unpack_dwrk(X,s,class)
%
% Extract the Butcher arrays of the solution from x::
%
%       x = [c' b' A ct' bt' At r] 
%
% (A,At are stored row-by-row)

[n,n2]=set_n_dwrk(s,class);

A=zeros(s); At=zeros(s);

switch class 
    case 'erk'
        c =[0 X(1:s-1)']';       b =X(s:2*s-1); 
        ct=[0 X(1+n2:s-1+n2)']'; bt=X(s+n2:2*s-1+n2); 
        for i=1:s
            A(i,1:i-1)=X(2*s+(i-2)*(i-1)/2   :2*s-1+i*(i-1)/2   );
            At(i,1:i-1)=X(2*s+(i-2)*(i-1)/2+n2:2*s-1+i*(i-1)/2+n2);
        end

    case 'dirk'
        c  = X(1:s);
        ct = X(1+n2:s+n2);
        b =X(s+1:2*s);
        bt=X(s+n2+1:2*s+n2);
        for i=1:s
            A(i,1:i)  = X(2*s+i*(i+1)/2-(i-1)   :2*s+i*(i+1)/2   );
            At(i,1:i) = X(2*s+i*(i+1)/2-(i-1)+n2:2*s+i*(i+1)/2+n2);
        end

    case 'irk'
        c =X(1:s);
        ct=X(1+n2:s+n2);
        b =X(s+1:2*s);
        bt=X(s+n2+1:2*s+n2);
        A =reshape(X(2*s+1:2*s+s^2),s,s)';
        At=reshape(X(n2+2*s+1:n2+2*s+s^2),s,s)';
end
end
