clearvars;clc;
% Parameters
pi = 0.8; nu_1 = 0.1; eta_1 = 0.2; mu_L = 0.05; theta_1=0.5; sigma_1=0.5;

R_d=(eta_1*pi)/((nu_1+mu_L)*(pi+mu_L));
R_s=R_d - (sigma_1^2*pi^2)/((4*theta_1)*(nu_1+mu_L)*(pi+mu_L)^2);

disp(['R_d = ',num2str(R_d)])
disp(['R_s = ',num2str(R_s)])