if diag_on
    diag.final_solve.solve = diag_solve; diag.final_solve.h = h; diag.final_solve.h_max = h_max;
    diag.final_solve.h_min = h_min; diag.final_solve.x_opt = x_opt; diag.final_solve.v = v;
    diag.final_solve.status = status;
    diag.n_iter = i_step;
end
