function n=set_n(s,class)
%function n=set_n(s,class)
% Set total number of decision variables
switch class
  case 'irk'  %Fully implicit
    n=2*s+s^2+1;
  case 'irk5'  %Implicit, p>=5 (so first row of A is zero)
    n=2*s+s*(s-1);
  case 'dirk' %Diagonally Implicit (A is lower triangular)
    n=2*s+s*(s+1)/2+1;
  case 'dirk5' %Diagonally Implicit, p>=5 (lower tri. and first row of A is zero)
    n=2*s+s*(s+1)/2-1;
  case 'sdirk' %Singly Diagonally Implicit
    n=2*s+s*(s+1)/2+1-s+1;
  case 'erk' % Explicit (A is strictly lower triangular)
    n=2*s+s*(s-1)/2;

  %====================
  % Low-storage classes
  %====================
  case '2S'             %Low-storage 2S methods of Ketcheson
    n = 3*s - 3;
  case '2Sstar'         %Low-storage 2S* methods of Ketcheson
    n = 3*s - 3;
  case '2Semb'          %Low-storage 2S* embedded pairs of Ketcheson
    n = 3*s - 1;
  case '3Sstar'         %Low-storage 3S* methods of Ketcheson
    n = 4*s - 6;
  case '3Sstaremb'      %Low-storage 3S* embedded pairs of Ketcheson
    n = 4*s - 3;
end
