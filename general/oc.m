function tau = oc(x,class,s,p)
% function tau=oc(x,class,s,p)
% Order conditions for RK methods

global oc_form

% Get Butcher coefficients
[A,b,c] = unpack_rk(x,s,class);

if strcmp(oc_form,'albrecht')
    tau = oc_albrecht(A,b,c,p);
elseif strcp(oc_form,'butcher')
    tau = oc_butcher(A,b,c,p);
end
