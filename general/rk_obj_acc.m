function [r]=rk_obj_acc(x,class,s,p)
%function [r,g]=rk_obj_acc(x,class,s,p)

r = errcoeff(x,class,s,p);