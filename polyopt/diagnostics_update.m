    if diag_on
        diag.step(i_step).solve = diag_solve; diag.step(i_step).h = h; diag.step(i_step).h_max = h_max;
        diag.step(i_step).h_min = h_min; diag.step(i_step).x_opt = x_opt; diag.step(i_step).v = v;
        diag.step(i_step).status = status;

        if strcmp(status,'Solved') || strcmp(status,'Inaccurate/Solved')
            diag.failed_steps(end+1) = i_step;
        end
 
    end


