function p = check_RK_order(A,b,c)

% Determines order of a RK method, up to sixth order

% For an s-stage method, input A should be a s x s matrix;
% b and c should are column vectors of length s

eps = 1.e-14;
s = length(b); % number of stages
unity = ones(s,1);
p = 0;

% order 1
t(1) = b'*unity - 1;
if max(abs(t))<eps
  p = 1;
end

% order 2
t(1) = b'*c - 1/2;
if (abs(t)<eps && p==1)
    p = 2;
end

% order 3
t(1) = b'*diag(c)*c - 1/3;
t(2) = b'*A*c - 1/6;
if(max(abs(t))<eps && p==2)
    p=3;
end

% order 4
t(1) = b'*diag(c)^2*c - 1/4;
t(2) = b'*diag(A*c)*c - 1/8;
t(3) = b'*A*diag(c)*c - 1/12;
t(4) = b'*A^2*c - 1/24;
if(max(abs(t))<eps && p==3)
    p=4;
end

% order 5
t(1) = b'*diag(c)^3*c - 1/5;
t(2) = b'*diag(c)^2*A*c - 1/10;
t(3) = b'*diag(c)*A*diag(c)*c - 1/15;
t(4) = b'*diag(c)*A^2*c - 1/30;
t(5) = b'*diag(A*c)*A*c - 1/20;
t(6) = b'*A*diag(c)^2*c - 1/20;
t(7) = b'*A*diag(c)*A*c - 1/40;
t(8) = b'*A^2*diag(c)*c - 1/60;
t(9) = b'*A^3*c - 1/120;
if(max(abs(t))<eps && p==4)
    p=5;
end


% order 6
t(1) = b'*diag(c)^4*c - 1/6;
t(2) = b'*diag(c)^3*A*c - 1/12;
t(3) = b'*diag(c)*diag(A*c)*A*c - 1/24;
t(4) = b'*diag(c)^2*A*diag(c)*c - 1/18;
t(5) = b'*diag(A*c)*A*diag(c)*c - 1/36;
t(6) = b'*diag(c)*A*diag(c)^2*c - 1/24;
t(7) = b'*A*diag(c)^3*c - 1/30;
t(8) = b'*diag(c)^2*A^2*c - 1/36;
t(9) = b'*diag(A*c)*A^2*c - 1/72;
t(10) = b'*diag(c)*A*diag(c)*A*c - 1/48;
t(11) = b'*A*diag(c)^2*A*c - 1/60;
t(12) = b'*A*diag(A*c)*A*c - 1/120;
t(13) = b'*diag(c)*A^2*diag(c)*c - 1/72;
t(14) = b'*A*diag(c)*A*diag(c)*c - 1/90;
t(15) = b'*A^2*diag(c)^2*c - 1/120;
t(16) = b'*diag(c)*A^3*c - 1/144;
t(17) = b'*A*diag(c)*A^2*c - 1/180;
t(18) = b'*A^2*diag(c)*A*c - 1/240;
t(19) = b'*A^3*diag(c)*c - 1/360;
t(20) = b'*A^4*c - 1/720;
if(max(abs(t))<eps && p==5)
    p=6;
    disp('This method has order at least 6');
end

end
