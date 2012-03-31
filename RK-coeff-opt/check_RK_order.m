function p = check_RK_order(A,b,c,problem_class)
%function p = check_RK_order(A,b,c)
% Determines order of a RK method, up to sixth order
% For an s-stage method, input A should be a s x s matrix;
% b and c should are column vectors of length s

done=0;
p=0;
eps = 1.e-14;

if strcmp(problem_class,'nonlinear')
     % order 1
    if abs(sum(b)-1)<eps
      p = 1;
    else
      done=1;
    end
    while done==0
        cond=oc_butcher(A,b,c,p+1);
        if max(abs(cond))<eps
            p=p+1;
        else
            done=1;
        end
    end
else
    while done==0
        linear_cond = linear_order_conditions(A,b,p+1);
        if abs(linear_cond)<eps
            p=p+1;
        else
            done=1;
        end
    end
end
        
    
