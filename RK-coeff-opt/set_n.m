function n=set_n(s,class,k)
% function n=set_n(s,class)
% Set total number of decision variables

switch class
    %=====================
    % RK classes
    %=====================
    case 'erk' % Explicit (A is strictly lower triangular)
      n=2*s+s*(s-1)/2;
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

    %=====================
    % Multistep RK classes
    %=====================
    case 'emsrk1' 
        n=.5*(s^2-s)+s*k+k;         %Explicit Type1
    case 'emsrk2' 
        n=.5*(s^2-s) +2*s*k-s +1; %Explicit Type2
    case 'imsrk1'  
        n=s^2+s*k+k;              %Fully implicit Type1
    case 'imsrk2'  
        n=s^2+2*s*k-s+1;          %Fully implicit Type2
    case 'dimsrk1'                       
        n=.5*(s^2+s)+s*k+k;         %Diagonally Implicit Type1
    case 'dimsrk2'                       
        n=.5*(s^2+s)+2*s*k-s+1      %Diagonally Implicit Type2
end
