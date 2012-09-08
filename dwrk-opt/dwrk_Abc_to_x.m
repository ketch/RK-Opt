function X=dwrk_Abc_to_x(A,At,b,bt,c,ct,r,class)
% Create vector X for optrk from A,b,c,r and class
% for downwind RK methods

X=[];
s=length(b)
switch class
  case 'erk'
    n=2*(2*s+s*(s-1)/2)-1; n2=(n-1)/2;
    for i=2:s
      X(2*s+(i-2)*(i-1)/2:2*s-1+i*(i-1)/2)=A(i,1:i-1);
      X(2*s+(i-2)*(i-1)/2+n2:2*s-1+i*(i-1)/2+n2)=At(i,1:i-1);
    end
    X(s:2*s-1)=b';
    X(s+n2:2*s-1+n2)=bt';
    X(1:s-1)=c(2:s)';
    X(1+n2:s-1+n2)=ct(2:s)';
    X(n)=-r;
end
