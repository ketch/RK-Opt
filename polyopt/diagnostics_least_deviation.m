if diag_on
    diag.b = b;
    diag.c = c;
    diag.fixed_coefficients = fixed_coefficients;
    diag.cond_b = cond(b);
    diag.cond_c = cond(c);
    diag.h = h;
    diag.precision = precision;
    diag.poly_coeffs = poly_coeffs;
    diag.cvx.cputime = cvx_cputime;
    diag.cvx.slvitr = cvx_slvitr;
    diag.cvx.slvtol = cvx_slvtol;
    diag.cvx.status = cvx_status;
else 
    diag = [];
end
