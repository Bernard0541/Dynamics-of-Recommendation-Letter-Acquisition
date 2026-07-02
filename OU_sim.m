clear; clc; close all;

%% Time settings
Tend = 100;
dt = 0.001;
t = 0:dt:Tend;
n = length(t);

%% Initial conditions
SU = zeros(1,n); SI = zeros(1,n); C = zeros(1,n); R = zeros(1,n); F = zeros(1,n);
T = zeros(1,n); A = zeros(1,n); U = zeros(1,n); N = zeros(1,n);

SU(1)=2000; SI(1)=1100; C(1)=50; R(1)=200; F(1)=40;
T(1)=1000; A(1)=30; U(1)=100; N(1)=90;

%% Parameters
LambdaS = 20; LambdaL = 10;
KS = 5000; KL = 2000;

rho = 0.4; kappa = 0.5;
alpha1 = 0.5; alpha2 = 0.1;

gamma1 = 0.1; gamma2 = 0.2; gamma3 = 0.3;
omega = 0.7; piL = 0.8;
epsilon = 0.5;

nu1 = 0.1; nu2 = 0.1;
muS = 0.005; muL = 0.05;

lambda = 1;   % set your value

%% Mean OU parameters
bar_beta1 = 0.8; beta10 = 0.8;
bar_beta2 = 0.5; beta20 = 0.5;
bar_beta3 = 0.7; beta30 = 0.7;

bar_eta1 = 0.2; eta10 = 0.2;
bar_eta2 = 0.3; eta20 = 0.3;

%% OU noise parameters
sigma_beta1 = 0.5;
sigma_beta2 = 0.5;
sigma_beta3 = 0.5;
sigma_eta1  = 0.5;
sigma_eta2  = 0.5;

theta_beta1 = 0.5;
theta_beta2 = 0.5;
theta_beta3 = 0.5;
theta_eta1  = 0.5;
theta_eta2  = 0.5;

%% Euler-Maruyama scheme
for i = 1:n-1
    ti = t(i);

    P = SU(i) + SI(i) + C(i) + R(i) + F(i);
    L = T(i) + A(i) + U(i) + N(i);

    b1 = bar_beta1 + (beta10 - bar_beta1)*exp(-theta_beta1*ti);
    b2 = bar_beta2 + (beta20 - bar_beta2)*exp(-theta_beta2*ti);
    b3 = bar_beta3 + (beta30 - bar_beta3)*exp(-theta_beta3*ti);

    e1 = bar_eta1 + (eta10 - bar_eta1)*exp(-theta_eta1*ti);
    e2 = bar_eta2 + (eta20 - bar_eta2)*exp(-theta_eta2*ti);

    sB1 = sigma_beta1/sqrt(2*theta_beta1)*sqrt(1 - exp(-2*theta_beta1*ti));
    sB2 = sigma_beta2/sqrt(2*theta_beta2)*sqrt(1 - exp(-2*theta_beta2*ti));
    sB3 = sigma_beta3/sqrt(2*theta_beta3)*sqrt(1 - exp(-2*theta_beta3*ti));

    sE1 = sigma_eta1/sqrt(2*theta_eta1)*sqrt(1 - exp(-2*theta_eta1*ti));
    sE2 = sigma_eta2/sqrt(2*theta_eta2)*sqrt(1 - exp(-2*theta_eta2*ti));

    dWb1 = sqrt(dt)*randn;
    dWb2 = sqrt(dt)*randn;
    dWb3 = sqrt(dt)*randn;
    dWe1 = sqrt(dt)*randn;
    dWe2 = sqrt(dt)*randn;

    SU(i+1) = SU(i) + ((1-rho)*LambdaS*(1-P/KS) ...
        - (alpha2*A(i) + gamma1*lambda + b2*A(i) + kappa + muS)*SU(i))*dt ...
        - sB2*A(i)*SU(i)*dWb2;

    SI(i+1) = SI(i) + (rho*LambdaS*(1-P/KS) + kappa*SU(i) ...
        - (alpha1*A(i) + b3*A(i) + gamma2*lambda + muS)*SI(i))*dt ...
        - sB3*A(i)*SI(i)*dWb3;

    C(i+1) = C(i) + ((alpha1*SI(i) + alpha2*SU(i))*A(i) ...
        - (gamma3*lambda + b1*A(i) + muS)*C(i))*dt ...
        - sB1*A(i)*C(i)*dWb1;

    R(i+1) = R(i) + ((b1*C(i) + b2*SU(i) + b3*SI(i))*A(i) ...
        - muS*R(i))*dt ...
        + sB1*C(i)*A(i)*dWb1 ...
        + sB2*SU(i)*A(i)*dWb2 ...
        + sB3*SI(i)*A(i)*dWb3;

    F(i+1) = F(i) + (gamma1*SU(i) + gamma2*SI(i) + gamma3*C(i))*lambda*dt;

    T(i+1) = T(i) + (LambdaL*(1-L/KL) - (piL + muL)*T(i))*dt;

    A(i+1) = A(i) + ((1-omega)*piL*T(i) + nu2*U(i) + nu1*N(i) ...
        - (epsilon + muL + e1*lambda)*A(i))*dt ...
        - sE1*lambda*A(i)*dWe1;

    U(i+1) = U(i) + (omega*piL*T(i) + epsilon*A(i) ...
        - (nu2 + muL + e2)*U(i))*dt ...
        - sE2*U(i)*dWe2;

    N(i+1) = N(i) + (e1*lambda*A(i) + e2*U(i) ...
        - (nu1 + muL)*N(i))*dt ...
        + sE1*lambda*A(i)*dWe1 ...
        + sE2*U(i)*dWe2;

    % Avoid negative populations
    SU(i+1)=max(SU(i+1),0); SI(i+1)=max(SI(i+1),0);
    C(i+1)=max(C(i+1),0); R(i+1)=max(R(i+1),0);
    F(i+1)=max(F(i+1),0); T(i+1)=max(T(i+1),0);
    A(i+1)=max(A(i+1),0); U(i+1)=max(U(i+1),0);
    N(i+1)=max(N(i+1),0);
end

%% Plot results
figure;
plot(t,SU,'LineWidth',1.5); hold on;
plot(t,SI,'LineWidth',1.5);
plot(t,C,'LineWidth',1.5);
plot(t,R,'LineWidth',1.5);
legend('S_U','S_I','C','R');
xlabel('Time'); ylabel('Population');
title('Human population classes');

figure;
plot(t,F,'LineWidth',1.5);
hold on;
legend('F');
xlabel('Time'); ylabel('Population');
title('Human population classes');

figure;
plot(t,T,'LineWidth',1.5); hold on;
plot(t,A,'LineWidth',1.5);
plot(t,U,'LineWidth',1.5);
plot(t,N,'LineWidth',1.5);
legend('T','A','U','N');
xlabel('Time'); ylabel('Population');
title('Livestock/vector population classes');