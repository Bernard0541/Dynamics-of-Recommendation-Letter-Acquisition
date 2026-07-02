clear; clc; close all;

%% Parameters
pi_val = 0.8;
muL = 0.05;

sigma_eta1 = 0.5;
theta_eta1 = 0.5;

%% Parameter ranges
eta1 = linspace(0,1,150);
nu1  = linspace(0,1,150);

[ETA1, NU1] = meshgrid(eta1, nu1);

r = pi_val/(pi_val + muL);

%% R0 surfaces
R0d = (ETA1 .* r) ./ (NU1 + muL);

R0s = (ETA1 .* r ...
    - (sigma_eta1^2/(4*theta_eta1))*r^2) ./ (NU1 + muL);

%% ===============================
%% FIGURE 1: Contour plots (RECOMMENDED)
%% ===============================
figure;

subplot(1,2,1)
contourf(ETA1, NU1, R0d, 20); hold on
contour(ETA1, NU1, R0d, [1 1], 'r', 'LineWidth', 3) % threshold line
colorbar
xlabel('\eta_1')
ylabel('\nu_1')
title('(a) \mathcal{R}_0^d(\eta_1,\nu_1)')
set(gca,'FontSize',12,'LineWidth',2)
grid on

subplot(1,2,2)
contourf(ETA1, NU1, R0s, 20); hold on
contour(ETA1, NU1, R0s, [1 1], 'r', 'LineWidth', 3)
colorbar
xlabel('\eta_1')
ylabel('\nu_1')
title('(b) \mathcal{R}_0^s(\eta_1,\nu_1)')
set(gca,'FontSize',12,'LineWidth',2)
grid on

figure;
surf(ETA1, NU1, R0d)
xlabel('\eta_1')
ylabel('\nu_1')
zlabel('\mathcal{R}_0^d')
title('Surface of R0^d')
shading interp
colorbar