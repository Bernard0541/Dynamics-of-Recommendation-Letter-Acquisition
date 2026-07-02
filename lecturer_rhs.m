function dYL = lecturer_rhs(~,YL,p)
    T = YL(1); A = YL(2); U = YL(3); I = YL(4);

    L = T + A + U + I;
    lambda = (U + I)/L;

    dT  = p.Lambda_L*(1 - L/p.K_L) - (p.pi + p.mu_L)*T;
    dA  = (1-p.omega)*p.pi*T + p.nu_2*U + p.nu_1*I - (p.eta_1*lambda + p.epsilon + p.mu_L)*A;
    dU  = p.omega*p.pi*T + p.epsilon*A - (p.nu_2 + p.eta_2 + p.mu_L)*U;
    dI  = p.eta_1*lambda*A + p.eta_2*U - (p.nu_1 + p.mu_L)*I;

    dYL = [dT; dA; dU; dI];
end