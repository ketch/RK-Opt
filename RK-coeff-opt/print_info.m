function info = print_info(status,p,order)
%info = display_info(status,p,order)
    if (status>0 && p==order)
        fprintf('\n===========================')
        fprintf('\nConverged to a solution. \n\n')
        fprintf('The RK coefficients satisfy the order conditions. \n')
        fprintf('Order of accuracy: %d \n', order)
        fprintf('===========================\n')
    
    else
        fprintf('\n===========================\n')
        if (status == 0)
            fprintf('Too many function evaluations or iterations. \n\n');
        elseif (status == -1)
            fprintf('Stopped by output/plot function. \n\n');
        elseif (status == -2)
            fprintf('No feasible point found. \n\n');
        elseif (status == -3)
            fprintf('Problem seems unbounded. \n\n');
        end
 
    fprintf('The RK method found does not satisfy the order conditions. \n');
    fprintf('===========================\n')
end
    
    
info = 1;
end

