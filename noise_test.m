clear; clc; close all;

%% Parameters
eta1 = 0.2;
pi_val = 0.8;
nu1 = 0.1;
muL = 0.05;
%theta_eta1 = 0.5;
sigma_eta1 = 0.5;

%% Range of noise intensity
theta_eta1 = linspace(0,1,100);

%% Common ratio
r = pi_val/(pi_val + muL);

%% Deterministic and stochastic indices
R0d = (eta1*r)/(nu1 + muL);

R0s = R0d - (sigma_eta1.^2 ./ (4*theta_eta1*(nu1 + muL))) .* r^2;

%% Plot
figure;
plot(theta_eta1, R0s, 'b-', 'LineWidth', 5)
hold on
yline(1, 'r--', 'LineWidth', 5);
yline(R0d, 'k:', 'LineWidth', 5);

xlabel('$\theta_{\eta_1}$')
ylabel('$\mathcal{R}_0^s$')
legend('$\mathcal{R}_0^s$', 'Threshold = 1', '$\mathcal{R}_0^d$', 'Location', 'best')
%title('Effect of noise intensity on stochastic institutional discouragement index')
grid on; box on;
set(gca, 'LineWidth', 2)