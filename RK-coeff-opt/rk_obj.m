function [r,g]=rk_obj(x,class,s,p,objective)
%function [r,g]=rk_obj(x,class,s,p,objective)

if strcmp(objective,'ssp')
    r=x(end);
    g=zeros(size(x));
    g(end)=1;
elseif strcmp(objective,'acc') % Will fail for multistep-RK methods
    [A,b,c]=unpack_rk(x,s,class);
    r=errcoeff(A,b,c,p);
    g=[];
end
