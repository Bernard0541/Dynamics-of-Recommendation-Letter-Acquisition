clearvars; clc;
Y0 = [2000, 1100, 50, 200, 40, 1000, 30, 100, 90]; 
p.Lambda_S = 20.0;
p.Lambda_L = 10.0;
p.K_S = 5000;
p.K_L = 2000;       
p.rho = 0.4;
p.kappa = 0.5;
p.alpha_1 = 0.5;
p.alpha_2 = 0.1;
p.gamma_1 = 0.1;
p.gamma_2 = 0.2;
p.gamma_3 = 0.3;
p.beta_1 = 0.8;
p.beta_2 = 0.5;
p.beta_3 = 0.7;
p.mu_S = 0.005;
p.mu_L = 0.05;
p.pi = 0.8;
p.omega = 0.7;
p.epsilon = 0.5;
p.nu_1 = 0.1;
p.nu_2 = 0.1;
p.eta_1 = 0.2;
p.eta_2 = 0.3;

tspan = [0 400]; 
[t, Y] = ode45(@(t, Y) recommendation_system(t, Y, p), tspan, Y0);

figure(1);
plot(t, Y(:,1), 'LineWidth', 5);
hold on
plot(t, Y(:,2), 'LineWidth', 5);
plot(t, Y(:,3), 'LineWidth', 5);
plot(t, Y(:,4), 'LineWidth', 5);
legend('Uninformed','Informed','Connected','Received','Location','best')
xlabel('Time (Years)')
ylabel('Number of Students')
box on; 
set(gca,'LineWidth',2);
xlim([0 2])
grid on;
hold off

figure(2);
plot(t, Y(:,6), 'LineWidth', 5);
hold on
plot(t, Y(:,7), 'LineWidth', 5);
plot(t, Y(:,8), 'LineWidth', 5);
plot(t, Y(:,9), 'LineWidth', 5);
legend('Recruited','Available','Unavailable','Negative','Location','best')
xlabel('Time (Years)')
ylabel('Number of Lecturers')
xlim([0 50])
box on;
set(gca,'LineWidth',2);
grid on;
hold off

figure(3);
plot(t, Y(:,5), 'LineWidth', 5);
hold on
%legend('F');
xlabel('Time (Years)')
ylabel('Cumulative Student failures')
xlim([0 50])
box on; 
set(gca,'LineWidth',2);
grid on;
hold off

function dY = recommendation_system(t, Y, p)
    % Assign state variables for readability
    SU=Y(1); SI=Y(2); C=Y(3); R=Y(4); F=Y(5); T=Y(6); A=Y(7); U=Y(8); N=Y(9);
    
    % Update Total Populations and Force of Infection dynamically
    P = SU + SI + C + R;
    L = T + A + U + N;
    lambda = (U + N) / L;
    
    % Equations
    dSU = (1-p.rho)*p.Lambda_S*(1 - P/p.K_S) - (p.alpha_2*A + p.gamma_1*lambda + p.beta_2*A + p.kappa + p.mu_S)*SU;
    dSI = p.rho*p.Lambda_S*(1 - P/p.K_S) - (p.alpha_1*A + p.beta_3*A + p.gamma_2*lambda + p.mu_S)*SI + p.kappa*SU;
    dC  = (p.alpha_1*SI + p.alpha_2*SU)*A - (p.gamma_3*lambda + p.beta_1*A + p.mu_S)*C;
    dR  = (p.beta_1*C + p.beta_2*SU + p.beta_3*SI)*A - p.mu_S*R;
    
    dF  = (p.gamma_1*SU + p.gamma_2*SI + p.gamma_3*C)*lambda;
    
    dT  = p.Lambda_L*(1 - L/p.K_L) - (p.pi + p.mu_L)*T;
    dA  = (1-p.omega)*p.pi*T + p.nu_2*U + p.nu_1*N - (p.eta_1*lambda + p.epsilon + p.mu_L)*A;
    dU  = p.omega*p.pi*T + p.epsilon*A - (p.nu_2 + p.eta_2 + p.mu_L)*U;
    dN  = p.eta_1*lambda*A + p.eta_2*U - (p.nu_1 + p.mu_L)*N;
    
    dY = [dSU; dSI; dC; dR; dF; dT; dA; dU; dN];
end