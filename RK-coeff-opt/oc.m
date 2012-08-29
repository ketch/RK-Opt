function tau=oc(x,class,s,p,Aeq,beq)
% function tau=oc(x,class,s,p)
% Order conditions for RKMs
% This is just a small wrapper

oc_form = 'albrecht';

[A,b,c]=unpack_rk(x,s,class);

if strcmp(oc_form,'albrecht')
    tau = oc_albrecht(A,b,c,p);
elseif strcp(oc_form,'butcher')
    tau = oc_butcher(A,b,c,p);
end

lin_tau = Aeq*x'-beq';
tau = [tau lin_tau'];
