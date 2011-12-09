function [r,g]=rk_obj(x,class,s,p,objective)
%function [r,g]=rk_obj(x,class,s,p,objective)

if strcmp(objective,'ssp')
    r=x(end);
    g=zeros(size(x));
    g(end)=1;
elseif strcmp(objective,'acc')
    r=errcoeff(x,class,s,p);
    g=[];
end
