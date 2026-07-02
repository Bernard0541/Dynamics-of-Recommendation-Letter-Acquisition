function forward_bif
clc; clearvars; close all;

%% Initial conditions
Y0 = [2000; 1100; 50; 200; 40; 1000; 30; 100; 90];
% [SU SI C R F T A U N]

%% Parameters
p.Lambda_S = 20.0;
p.Lambda_L = 10.0;
p.K_S = 5000;
p.K_L = 2000;
p.rho = 0.4;
p.kappa = 0.5;
p.alpha_1 = 0.5;
p.alpha_2 = 0.1;
p.beta_1 = 0.8;
p.beta_2 = 0.5;
p.beta_3 = 0.7;
p.gamma_1 = 0.1;
p.gamma_2 = 0.2;
p.gamma_3 = 0.3;
p.pi = 0.8;
p.eta_2 = 0.3;
p.nu_1 = 0.1;
p.nu_2 = 0.1;
p.mu_S = 0.005;
p.mu_L = 0.05;

%% Important for forward bifurcation about the SLE
p.omega = 0;
p.epsilon = 0;

%% Critical eta_1 value
eta1_star = ((p.nu_1 + p.mu_L)*(p.pi + p.mu_L))/p.pi;
fprintf('Critical eta_1 = %.6f\n', eta1_star);

%% Sweep eta_1
eta1_vals = linspace(0,1,300);
Nstar = zeros(size(eta1_vals));
Ustar = zeros(size(eta1_vals));
Astar = zeros(size(eta1_vals));
R0d_vals = zeros(size(eta1_vals));

tspan = [0 1000];

for j = 1:length(eta1_vals)

    p.eta_1 = eta1_vals(j);

    R0d_vals(j) = (p.eta_1*p.pi)/((p.nu_1+p.mu_L)*(p.pi+p.mu_L));

    [~,Y] = ode45(@(t,Y) recommendation_model(t,Y,p), tspan, Y0);

    Yend = Y(end,:);

    Ustar(j) = Yend(8);
    Nstar(j) = Yend(9);
    Astar(j) = Yend(7);
end

%% Figure 1: Forward bifurcation diagram
% figure;
% plot(R0d_vals, Nstar, 'b-', 'LineWidth', 3);
% hold on
% plot(R0d_vals, zeros(size(R0d_vals)), 'k--', 'LineWidth', 2)
% xline(1, 'r--', 'LineWidth', 2)
% 
% xlabel('\mathcal{R}_0^d')
% ylabel('Equilibrium negative influencers, N^*')
% title('Forward bifurcation diagram')
% legend('Unsupportive lecturer equilibrium branch', 'Supportive lecturer equilibrium branch',...
%     '\mathcal{R}_0^d=1 threshold', 'Location','best')
% grid on; box on;
% set(gca,'FontSize',12,'LineWidth',2)

%% Figure 1: Forward bifurcation diagram

figure;
hold on

% Split branches by threshold R0d = 1
idx_stable_SLE   = R0d_vals < 1;
idx_unstable_SLE = R0d_vals > 1;
idx_stable_ULE   = R0d_vals > 1;

% Stable SLE: blue line, N*=0 for R0d < 1
plot(R0d_vals(idx_stable_SLE), zeros(size(R0d_vals(idx_stable_SLE))), ...
    'b-', 'LineWidth', 5);

% Unstable SLE: green dashed line, N*=0 for R0d > 1
plot(R0d_vals(idx_unstable_SLE), zeros(size(R0d_vals(idx_unstable_SLE))), ...
    'g-', 'LineWidth', 5);

% Stable ULE: red line, N*>0 for R0d > 1
plot(R0d_vals(idx_stable_ULE), Nstar(idx_stable_ULE), ...
    'r-', 'LineWidth', 5);

% Threshold line
xline(1, 'k:', 'LineWidth', 4);

xlabel('$\mathcal{R}_0^d$')
ylabel('Equilibrium negative-influencers $(N^*)$')
%title('Forward bifurcation diagram')

legend('Stable $\mathcal{E}^0$','Unstable $\mathcal{E}^0$','Stable $\mathcal{E}^*$','$\mathcal{R}_0^d=1$','Location','best')
xlim([0 6])
grid on; box on;
set(gca,'LineWidth',2)
hold off


end

%% Model
function dY = recommendation_model(~,Y,p)

SU = Y(1);
SI = Y(2);
C  = Y(3);
R  = Y(4);
F  = Y(5);
T  = Y(6);
A  = Y(7);
U  = Y(8);
N  = Y(9);

P = SU + SI + C + R;
L = T + A + U + N;
lambda = (U + N)/max(L,1e-12);

dY = zeros(9,1);

dY(1) = (1-p.rho)*p.Lambda_S*(1 - P/p.K_S) ...
      - (p.alpha_2*A + p.gamma_1*lambda + p.beta_2*A + p.kappa + p.mu_S)*SU;

dY(2) = p.rho*p.Lambda_S*(1 - P/p.K_S) ...
      - (p.alpha_1*A + p.beta_3*A + p.gamma_2*lambda + p.mu_S)*SI ...
      + p.kappa*SU;

dY(3) = (p.alpha_1*SI + p.alpha_2*SU)*A ...
      - (p.gamma_3*lambda + p.beta_1*A + p.mu_S)*C;

dY(4) = (p.beta_1*C + p.beta_2*SU + p.beta_3*SI)*A ...
      - p.mu_S*R;

dY(5) = (p.gamma_1*SU + p.gamma_2*SI + p.gamma_3*C)*lambda;

dY(6) = p.Lambda_L*(1 - L/p.K_L) ...
      - (p.pi + p.mu_L)*T;

dY(7) = (1-p.omega)*p.pi*T + p.nu_2*U + p.nu_1*N ...
      - (p.eta_1*lambda + p.epsilon + p.mu_L)*A;

dY(8) = p.omega*p.pi*T + p.epsilon*A ...
      - (p.nu_2 + p.eta_2 + p.mu_L)*U;

dY(9) = p.eta_1*lambda*A + p.eta_2*U ...
      - (p.nu_1 + p.mu_L)*N;

end