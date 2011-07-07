function n=set_n(s,class,regs)
%function n=set_n(s,class,regs)
% Set total number of decision variables
switch regs
  case 2
    switch class
      case 'plain'
        n = 3*s - 3;
      case 'star'
        n = 3*s - 3;
      case 'emb'
        n = 3*s - 1;
    end

  case 3
    switch class
      case 'star'
        n = 4*s - 6;
      case 'staremb'
        n = 4*s - 3;
    end
end
