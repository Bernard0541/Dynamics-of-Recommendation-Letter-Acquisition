clear; clc; close all;

%% Parameters
%pi = 0.8; 
nu_1 = 0.1; 
eta_1 = 0.2;
mu_L = 0.05;
theta_1=0.5;
%sigma_1=0.5;

% Define ranges for eta1 and nu1
%eta_1_range = linspace(0.0, 1, 200); 
%nu_1_range = linspace(0.0, 1, 200);
pi_range = linspace(0.0, 1, 200);
%mu_L_range = linspace(0.0, 1, 200); 
%theta_1_range = linspace(0.0, 1, 200);
sigma_1_range = linspace(0.0, 1, 200);

% Initialize matrix for Rd and Rs
Rd = zeros(length(pi_range), length(sigma_1_range));
Rs = zeros(length(pi_range), length(sigma_1_range));

% Calculate R0 for each combination
for i = 1:length(pi_range)
    for j = 1:length(sigma_1_range)
        pi = pi_range(i);
        sigma_1 = sigma_1_range(j);

        Rd(i, j) = (eta_1*pi)/((nu_1+mu_L)*(pi+mu_L));
        Rs(i, j) = (eta_1*pi)/((nu_1+mu_L)*(pi+mu_L)) - (sigma_1^2*pi^2)/((4*theta_1)*(nu_1+mu_L)*(pi+mu_L)^2);
    end
end

figure(1);

% subplot(1,2,1)
% contourf(pi_range, sigma_1_range, Rd);
% hold on
% contour(pi_range, sigma_1_range, Rd, [1 1], 'r', 'LineWidth', 4) % threshold line
% colorbar
% xlabel('$\pi$')
% ylabel('$\sigma_1$')
% title('(a) $\mathcal{R}_0^d(\pi,\sigma_1)$')
% set(gca, 'LineWidth',2)
% grid on

%subplot(1,2,2)
contourf(pi_range, sigma_1_range, Rs);
hold on
contour(pi_range, sigma_1_range, Rs, [1 1], 'r', 'LineWidth', 4)
colorbar
xlabel('$\pi$')
ylabel('$\sigma_1$')
% title('(b) $\mathcal{R}_0^s(\pi,\sigma_1)$')
title('$\mathcal{R}_0^s(\pi,\sigma_1)$')
set(gca, 'LineWidth',2)
grid on