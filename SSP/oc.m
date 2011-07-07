function tau=oc(x,class,s,p)
%function tau=oc(x,class,s,p)
%Order conditions for RKMs
%This is just a small wrapper

[A,b,c]=unpack_rk(x,s,class);

tau = oc_albrecht(A,b,c);
