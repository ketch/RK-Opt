function [n,n2]=set_n_dwrk(s,class)
% function [n,n2]=set_n_dwrk(s,class)

switch class
  case 'erk' % Explicit (A, At are strictly lower triangular)
    n=2*(2*s+s*(s-1)/2)-1; % = s^2+3*s-1
  case 'dirk' % Diagonally implicit (A, At are lower triangular)
    n=2*(2*s+s*(s+1)/2)+1;
  case 'irk' % Implicit (A, At are full)
    n=2*s^2+4*s+1;
end

n2=(n-1)/2; %Storage required for just b,c,A
