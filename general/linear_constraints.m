function [Aeq,beq,lb,ub] = linear_constraints(s,class,n)
% function [Aeq,beq,lb,ub] = linear_constraints(s,class,n)
%
% Sets up:
% - the linear constraints, i.e. Aeq*x = beq
% - the upper and lower bounds on the unknowns, i.e. ub, lb

switch class
  case 'irk'       % Fully implicit
    Aeq = zeros(s+1,n); 
    beq = zeros(1,s+1);
    for i = 1:s    % Require c's to be row sums of A
      Aeq(i,i) = 1;
      Aeq(i,((2+(i-1))*s+1):((2+i)*s)) = -1;
    end
    Aeq(end,s+1:2*s) = 1; 
    beq(end) = 1;  % Require b's to sum to one
    lb = zeros(1,n); 
    lb(end) = -30;
    ub = 2+zeros(1,n); 
    ub(end) = 0; 
    ub(s+1:2*s) = 1;

  case 'irk5'    % Implicit, p>=5 (so first row of A is zero)
    Aeq = zeros(s,n); 
    beq = zeros(1,s);
    for i = 2:s  % Require c's to be row sums of A
      Aeq(i-1,i-1) = 1;
      Aeq(i-1,((2+(i-2))*s):((2+(i-1))*s-1)) = -1;
    end
    Aeq(end,s:2*s-1) = 1; 
    beq(end) = 1; % Require b's to sum to one
    lb = zeros(1,n); 
    lb(end) = -30;
    ub = 2+zeros(1,n); 
    ub(end) = 0; 
    ub(s:2*s-1) = 1;

  case 'dirk' % Diagonally Implicit (A is lower triangular)
    Aeq = zeros(s+1,n);
    beq = zeros(1,s+1);
    for i = 1:s
      Aeq(i,i) = 1;
      Aeq(i,2*s+1+i*(i-1)/2:2*s+i*(i+1)/2) = -1;
    end
    Aeq(end,s+1:2*s) = 1; 
    beq(end) = 1;
    lb = zeros(1,n); 
    lb(end) = -30;
    ub = 2+zeros(1,n); 
    ub(end) = 0; 
    ub(s+1:2*s) = 1;

  case 'dirk5' % Diagonally Implicit, p>=5 (lower tri. and first row of A is zero)
    Aeq = zeros(s,n);
    beq = zeros(1,s);
    for i = 2:s
      Aeq(i-1,i-1) = 1;
      Aeq(i-1,2*s-1+i*(i-1)/2:2*s-2+i*(i+1)/2) = -1;
    end
    Aeq(end,s:2*s-1) = 1; 
    beq(end) = 1;
    lb = zeros(1,n); 
    lb(end) = -30;
    ub = 2+zeros(1,n); 
    ub(end) = 0; 
    ub(s:2*s-1) = 1;

  case 'sdirk' % Singly Diagonally Implicit
    Aeq = zeros(s+1,n);
    beq = zeros(1,s+1);
    Aeq(1,1) = 1; 
    Aeq(1,2*s+1) = -1;
    for i = 2:s
      Aeq(i,i) = 1;
      Aeq(i,2*s+1) = -1;
      Aeq(i,2*s+2+(i-1)*(i-2)/2:2*s+1+i*(i-1)/2) = -1;
    end
    Aeq(end,s+1:2*s) = 1; 
    beq(end) = 1;
    lb = zeros(1,n); 
    lb(end) = -30;
    ub = 2+zeros(1,n); 
    ub(end) = 0; 
    ub(s+1:2*s) = 1;

  case 'erk' % Explicit (A is strictly lower triangular)
    Aeq = zeros(s,n);
    beq = zeros(1,s);
    for i = 2:s
      Aeq(i-1,i-1) = 1;
      Aeq(i-1,2*s+(i-2)*(i-1)/2:2*s-1+i*(i-1)/2) = -1;
    end
    Aeq(end,s:2*s-1) = 1; 
    beq(end) =1;
    lb = -5+zeros(1,n); 
    lb(end) = -30;
    ub = 5+zeros(1,n); 
    ub(end) = 0;% ub(s:2*s-1)=1;

  otherwise
    % for low-storage classes:
    Aeq = []; 
    beq =[];
    lb = -8+zeros(1,n); %lb(end)=-30;
    ub = 8+zeros(1,n); %ub(end)=0;

end

