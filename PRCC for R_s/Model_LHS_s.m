tic
clearvars;
close all;
clc

R_0 = 0.51672;

runs=1000000;
           
Parameter_settings_LHS;
pi_LHS=LHS_Call(0.8*0.9,pi,0.8*1.1,0,runs,'unif');
nu_1_LHS=LHS_Call(0.1*0.9,nu_1,0.1*1.1,0,runs,'unif');
eta_1_LHS=LHS_Call(0.2*0.9,eta_1,0.2*1.1,0,runs,'unif');
mu_L_LHS=LHS_Call(0.05*0.9,mu_L,0.05*1.1,0,runs,'unif');
sigma_1_LHS=LHS_Call(0.5*0.9,sigma_1,0.5*1.1,0,runs,'unif');
theta_1_LHS=LHS_Call(0.5*0.9,theta_1,0.5*1.1,0,runs,'unif');


LHSmatrix=[pi_LHS, nu_1_LHS, eta_1_LHS, mu_L_LHS, sigma_1_LHS, theta_1_LHS];


for x=1:runs 
    x;
    LHSmatrix (x,:);
    Rs = RRR(LHSmatrix,x);
    Rs_lhs(:,x)= [Rs];
end

E_R = mean(Rs_lhs);
Var_R = var(Rs_lhs);
disp(['E[R] = ', num2str(E_R)]);
disp(['Var[R] = ', num2str(Var_R)]);

figure(2);
histogram(Rs_lhs, 100, 'Normalization', 'pdf');
hold on;
xlabel('$\mathcal{R}_0^s$');
ylabel('Probability Density');
xline(E_R, 'Color', 'r', 'LineWidth', 2)
xline(R_0, 'Color', 'k', 'LineWidth', 2)
%xline(Var_R, 'Color', 'g', 'LineWidth', 2)
xline(prctile(Rs_lhs, 5), 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2)%, 'Label', 'R_0(5th)');
xline(prctile(Rs_lhs, 95), 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2)%, 'Label', 'R_0(95th)')
set(gca, 'LineWidth', 2);

%save Model_LHS.mat;

  [prcc sign sign_label]=PRCC(LHSmatrix,Rs_lhs,1,PRCC_var,0.05);

toc