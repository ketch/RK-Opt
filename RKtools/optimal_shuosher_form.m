function [v,alpha,beta] = optimal_shuosher_form(A,b,c)
%function [v,alpha,beta] = optimal_shuosher_form(A,b,c)


s = length(b);
r=am_radius(A,b,c);
K = [A;b'];
K(:,end+1)=0;
e = ones(s+1,1);

beta = K/(eye(s+1)+r*K);
alpha = r*beta;
v = (eye(s+1)-alpha)*e;

