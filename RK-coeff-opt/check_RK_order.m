function p = check_RK_order(A,b,c,problem_class)
% function p = check_RK_order(A,b,c)
% Determines order of a RK method, up to sixth order.
%
% For an s-stage method, input `A` should be a `s \times s` matrix;
% `b` and `c` should be column vectors of length `s`.

if nargin<4
    problem_class='nonlinear';
end

eps = 1.e-14;

if abs(sum(b)-1)<eps
    p = 1;
else
    p = 0;
    return
end

while p<10
    if strcmp(problem_class,'nonlinear')
        cond=oc_butcher(A,b,c,p+1);
    else
        cond = b'*A^(p)*ones(length(b),1) - 1/factorial(p+1);
    end
    if max(abs(cond))<eps
        p=p+1;
    else
        return
    end
end
