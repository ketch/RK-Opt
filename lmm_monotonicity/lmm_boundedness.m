function gamma = lmm_boundedness(a,b,N)
% Determine the largest step size under which the given linear multistep method
% is bounded.
% Based on work of Willem Hundsdorfer and collaborators.

gamma_max = 2.;
gamma_min = 0.;
prec = 0.01;

k = length(a);
rho = zeros(1,N+1);
ppi = zeros(1,N+1);

while gamma_max - gamma_min > prec

    gamma = 0.5 * (gamma_max + gamma_min);

    rho(1) = 1./(1.+gamma*b(1));
    ppi(1) = gamma*b(1)*rho(1);

    % Startup
    for n = 1:k
        rho(n+1) = (sum(a(1:n).*rho(n:-1:1)) - gamma*sum(b(2:n+1).*rho(n:-1:1)))/(1.+gamma*b(1));
    end

    % Main
    for n = k+1:N
        rho(n+1) = (sum(a.*rho(n:-1:n+1-k)) - gamma*sum(b(2:end).*rho(n:-1:n+1-k)) )/(1.+gamma*b(1));
        ppi(n+1) = gamma * sum(b.*rho(n+1:-1:n+1-k));
    end

    if all(rho>=0)
        gamma_min = gamma;
    else
        gamma_max = gamma;
    end


end
