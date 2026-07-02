function OU_parameter_sweep
clc; clearvars; close all;
rng(1);

% Time settings
Tend = 100;
dt = 0.001;
t = 0:dt:Tend;
n = length(t);

% Initial conditions
Y0 = [2000; 1100; 50; 200; 40; 1000; 30; 100; 90];
% [SU SI C R F T A U N]

% Fixed parameters
p.LambdaS = 20;
p.LambdaL = 10;
p.KS = 5000;
p.KL = 2000;
p.rho = 0.4;
p.kappa = 0.5;
p.alpha1 = 0.5;
p.alpha2 = 0.1;
p.gamma1 = 0.1;
p.gamma2 = 0.2;
p.gamma3 = 0.3;
p.omega = 0.7;
p.piL = 0.8;
p.epsilon = 0.5;
p.nu1 = 0.1;
p.nu2 = 0.1;
p.muS = 0.005;
p.muL = 0.05;

% OU mean values and initial values
p.bar_beta1 = 0.8; p.beta10 = 0.8;
p.bar_beta2 = 0.5; p.beta20 = 0.5;
p.bar_beta3 = 0.7; p.beta30 = 0.7;

p.bar_eta1 = 0.2; p.eta10 = 0.2; p.bar_eta2 = 0.3; p.eta20 = 0.3;

% OU noise parameters
p.sigma_beta1 = 0.5;
p.sigma_beta2 = 0.5;
p.sigma_beta3 = 0.5;
p.sigma_eta1  = 0.5;
p.sigma_eta2  = 0.5;

p.theta_beta1 = 0.5;
p.theta_beta2 = 0.5;
p.theta_beta3 = 0.5;
%p.theta_eta1  = 0.5;
p.theta_eta2  = 0.5;

% Parameter values to sweep
%piL_values = [0.10, 0.25, 0.50, 0.80, 1.00];
%eta_1_values = [0.20, 0.40, 0.60, 0.80, 1.00];
%nu_1_values = [0.10, 0.25, 0.50, 0.80, 1.00];
%mu_L_values = [0.01, 0.05, 0.10, 0.25, 0.50];
%sigma_eta_1_values = [0.10, 0.25, 0.50, 0.80, 1.00];
theta_eta_1_values = [0.10, 0.25, 0.50, 0.80, 1.00];

% Plot recovered students R(t)
figure(1);
hold on;

for j = 1:length(theta_eta_1_values)

    p.theta_eta1 = theta_eta_1_values(j);

    Y = simulate_OU_model(Y0,t,dt,p);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    R = Y(4,:);

    plot(t,R,'LineWidth',4, 'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])
end

xlabel('Time')
ylabel('Students who received recommendation letters ($R(t)$)')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

% Plot cumulative failures F(t)
figure(2);
hold on;

for j = 1:length(theta_eta_1_values)

    p.theta_eta1 = theta_eta_1_values(j);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    Y = simulate_OU_model(Y0,t,dt,p);

    F = Y(5,:);

    plot(t,F,'LineWidth',4,'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])

end

xlabel('Time')
ylabel('Cumulative failures ($F(t)$)')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

% Plot available and supportive lecturers A(t)
figure(3);
hold on;

for j = 1:length(theta_eta_1_values)
    p.theta_eta1 = theta_eta_1_values(j);
    Y = simulate_OU_model(Y0,t,dt,p);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    A = Y(7,:);

    plot(t,A,'LineWidth',4, 'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])

end

xlabel('Time')
ylabel('Available and supportive lecturers ($A(t)$)')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

% Plot discouragement population U(t) and N(t)
figure(4);
hold on;

for j = 1:length(theta_eta_1_values)
    p.theta_eta1 = theta_eta_1_values(j);
    Y = simulate_OU_model(Y0,t,dt,p);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    U = Y(8,:);
    N = Y(9,:);
    D = U + N;

    plot(t,D,'LineWidth',4, 'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])

end

xlabel('Time')
ylabel('Discouragement population ($U(t) + N(t)$)')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

% Plot Recruited or newly employed lecturers T(t)
figure(5);
hold on;

for j = 1:length(theta_eta_1_values)
    p.theta_eta1 = theta_eta_1_values(j);
    Y = simulate_OU_model(Y0,t,dt,p);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    T = Y(6,:);

    plot(t,T,'LineWidth',4, 'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])

end

xlabel('Time')
ylabel('Recruited lecturers ($T(t)$)')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

% Plot Uninformed students seeking recommendation letters SU(t)
figure(6);
hold on;

for j = 1:length(theta_eta_1_values)
    p.theta_eta1 = theta_eta_1_values(j);
    Y = simulate_OU_model(Y0,t,dt,p);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    SU = Y(1,:);

    plot(t,SU,'LineWidth',4, 'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])

end

xlabel('Time')
ylabel('Uninformed students ($S_U(t)$')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

% Plot Informed students seeking recommendation letters SI(t)
figure(7);
hold on;

for j = 1:length(theta_eta_1_values)
    p.theta_eta1 = theta_eta_1_values(j);
    Y = simulate_OU_model(Y0,t,dt,p);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    SI = Y(2,:);

    plot(t,SI,'LineWidth',4, 'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])

end

xlabel('Time')
ylabel('Informed students ($S_I(t)$)')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

% Plot connected students C(t)
figure(8);
hold on;

for j = 1:length(theta_eta_1_values)
    p.theta_eta1 = theta_eta_1_values(j);
    Y = simulate_OU_model(Y0,t,dt,p);

    R0d = (p.bar_eta1*p.piL)/((p.nu1+p.muL)*(p.piL+p.muL));
    R0s = R0d - (p.piL^2*p.sigma_eta1^2)/(4*p.theta_eta1*(p.nu1+p.muL)*(p.piL+p.muL)^2);

    C = Y(3,:);

    plot(t,C,'LineWidth',4, 'DisplayName',['$\theta_{\eta_1}$ = ' num2str(p.theta_eta1,'%.2f') ...
        ',  $\mathcal{R}_0^d$ = ' num2str(R0d,'%.3f') ',  $\mathcal{R}_0^s$ = ' num2str(R0s,'%.3f')])

end

xlabel('Time')
ylabel('Connected students ($C(t)$)')
legend('show','Location','best')
grid on; box on;
set(gca,'LineWidth',2)
hold off;

end
% OU stochastic model
function Y = simulate_OU_model(Y0,t,dt,p)
n = length(t); Y = zeros(9,n); Y(:,1) = Y0;
for i = 1:n-1
    ti = t(i);
    SU = Y(1,i); SI = Y(2,i); C  = Y(3,i); R  = Y(4,i); F  = Y(5,i); T  = Y(6,i); A  = Y(7,i); U  = Y(8,i); N  = Y(9,i);
    P = SU + SI + C + R; L = T + A + U + N;
    lambda = (U + N)/max(L,1e-12);

    b1 = p.bar_beta1 + (p.beta10 - p.bar_beta1)*exp(-p.theta_beta1*ti);
    b2 = p.bar_beta2 + (p.beta20 - p.bar_beta2)*exp(-p.theta_beta2*ti);
    b3 = p.bar_beta3 + (p.beta30 - p.bar_beta3)*exp(-p.theta_beta3*ti);

    e1 = p.bar_eta1 + (p.eta10 - p.bar_eta1)*exp(-p.theta_eta1*ti);
    e2 = p.bar_eta2 + (p.eta20 - p.bar_eta2)*exp(-p.theta_eta2*ti);

    sB1 = p.sigma_beta1/sqrt(2*p.theta_beta1)*sqrt(1-exp(-2*p.theta_beta1*ti));
    sB2 = p.sigma_beta2/sqrt(2*p.theta_beta2)*sqrt(1-exp(-2*p.theta_beta2*ti));
    sB3 = p.sigma_beta3/sqrt(2*p.theta_beta3)*sqrt(1-exp(-2*p.theta_beta3*ti));

    sE1 = p.sigma_eta1/sqrt(2*p.theta_eta1)*sqrt(1-exp(-2*p.theta_eta1*ti));
    sE2 = p.sigma_eta2/sqrt(2*p.theta_eta2)*sqrt(1-exp(-2*p.theta_eta2*ti));

    dWb1 = sqrt(dt)*randn; dWb2 = sqrt(dt)*randn; dWb3 = sqrt(dt)*randn; dWe1 = sqrt(dt)*randn; dWe2 = sqrt(dt)*randn;

    dSU = ((1-p.rho)*p.LambdaS*(1-P/p.KS) - (p.alpha2*A + p.gamma1*lambda + b2*A + p.kappa + p.muS)*SU)*dt - sB2*A*SU*dWb2;
    dSI = (p.rho*p.LambdaS*(1-P/p.KS) + p.kappa*SU - (p.alpha1*A + b3*A + p.gamma2*lambda + p.muS)*SI)*dt - sB3*A*SI*dWb3;
    dC = ((p.alpha1*SI + p.alpha2*SU)*A - (p.gamma3*lambda + b1*A + p.muS)*C)*dt - sB1*A*C*dWb1;
    dR = ((b1*C + b2*SU + b3*SI)*A - p.muS*R)*dt + sB1*C*A*dWb1 + sB2*SU*A*dWb2 + sB3*SI*A*dWb3;
    dF = (p.gamma1*SU + p.gamma2*SI + p.gamma3*C)*lambda*dt;
    dT = (p.LambdaL*(1-L/p.KL) - (p.piL + p.muL)*T)*dt;
    dA = ((1-p.omega)*p.piL*T + p.nu2*U + p.nu1*N - (p.epsilon + p.muL + e1*lambda)*A)*dt - sE1*lambda*A*dWe1;
    dU = (p.omega*p.piL*T + p.epsilon*A - (p.nu2 + p.muL + e2)*U)*dt - sE2*U*dWe2;
    dN = (e1*lambda*A + e2*U - (p.nu1 + p.muL)*N)*dt + sE1*lambda*A*dWe1 + sE2*U*dWe2;

    Y(:,i+1) = Y(:,i) + [dSU; dSI; dC; dR; dF; dT; dA; dU; dN];
    Y(:,i+1) = max(Y(:,i+1),0);
end
end