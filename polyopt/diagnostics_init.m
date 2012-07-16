if nargout >3           %Diagnostics requested
    if isempty(h_true)
        error('Enabling diagnostics assumes h_true is given!');
    end
    diag_on=true;
    diag.lam = lam; diag.s = s; diag.p = p; diag.basis = basis; diag.tol_bisect = tol_bisect;
    diag.tol_feasible = tol_feasible; diag.h_min = h_min; diag.h_max = h_max; diag.max_steps = max_steps; diag.h_true = h_true;
    diag.solved_steps = []; diag.failed_steps = []; diag.right_steps = []; diag.wrong_steps = [];
else
    diag_on=false;
end


