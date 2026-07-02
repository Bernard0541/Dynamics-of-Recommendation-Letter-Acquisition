function discouragement_survival
clc; clearvars; close all;
rng(1);

%% Initial conditions: [SU SI C R F T A U N]
Y0 = [2000; 1100; 50; 200; 40; 1000; 30; 100; 90];

%% Parameters
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

p.pi = 0.8; 
p.nu_1 = 0.1; 
p.nu_2 = 0.1;

p.mu_S = 0.005; 
p.mu_L = 0.05;

%% Supportive baseline
p.omega = 0;
p.epsilon = 0;

%% OU mean and initial values
p.beta1_bar = 0.8; 
p.beta2_bar = 0.5; 
p.beta3_bar = 0.7; 
p.eta2_bar  = 0.3;

p.beta10 = 0.8; 
p.beta20 = 0.5; 
p.beta30 = 0.7; 
p.eta20  = 0.3;

%% OU stochastic parameters
p.theta_beta1 = 0.5; 
p.theta_beta2 = 0.5; 
p.theta_beta3 = 0.5;
p.theta_eta1  = 0.5; 
p.theta_eta2  = 0.5;

p.sigma_beta1 = 0.10; 
p.sigma_beta2 = 0.10; 
p.sigma_beta3 = 0.10;
p.sigma_eta1  = 0.50; 
p.sigma_eta2  = 0.10;

%% Fast simulation settings
Tfinal = 100;
dt = 0.05;
nPaths = 10000;
Ntol = 1e-3;

%% LHS sample for eta_1
nEta = 500;
eta1_min = 0;
eta1_max = 1;

eta1_vals = lhsdesign(nEta,1);
eta1_vals = eta1_min + (eta1_max - eta1_min)*eta1_vals;
eta1_vals = sort(eta1_vals(:))';

%% Time grid
tgrid = 0:dt:Tfinal;
nSteps = length(tgrid);
sqrt_dt = sqrt(dt);

%% Storage
Psurv = zeros(size(eta1_vals));
R0s_vals = zeros(size(eta1_vals));
Nmean = zeros(size(eta1_vals));
Dmean = zeros(size(eta1_vals));   % D = U + N

%% Precompute diffusion amplitudes
phi_beta1 = p.sigma_beta1/sqrt(2*p.theta_beta1).*sqrt(1-exp(-2*p.theta_beta1*tgrid));
phi_beta2 = p.sigma_beta2/sqrt(2*p.theta_beta2).*sqrt(1-exp(-2*p.theta_beta2*tgrid));
phi_beta3 = p.sigma_beta3/sqrt(2*p.theta_beta3).*sqrt(1-exp(-2*p.theta_beta3*tgrid));
phi_eta1  = p.sigma_eta1 /sqrt(2*p.theta_eta1 ).*sqrt(1-exp(-2*p.theta_eta1*tgrid));
phi_eta2  = p.sigma_eta2 /sqrt(2*p.theta_eta2 ).*sqrt(1-exp(-2*p.theta_eta2*tgrid));

%% Main loop over eta_1
for j = 1:length(eta1_vals)

    eta1_bar = eta1_vals(j);
    eta10 = eta1_bar;

    R0d = (eta1_bar*p.pi)/((p.nu_1+p.mu_L)*(p.pi+p.mu_L));

    R0s_vals(j) = R0d ...
        - (p.sigma_eta1^2/(4*p.theta_eta1*(p.nu_1+p.mu_L))) ...
        * (p.pi/(p.pi+p.mu_L))^2;

    %% Simulate all paths at once
    Y = repmat(Y0,1,nPaths);

    for k = 1:nSteps

        t = tgrid(k);

        SU = Y(1,:); 
        SI = Y(2,:); 
        C  = Y(3,:); 
        R  = Y(4,:);
        T  = Y(6,:); 
        A  = Y(7,:); 
        U  = Y(8,:); 
        N  = Y(9,:);

        Ppop = SU + SI + C + R;
        Lpop = T + A + U + N;

        lambda = (U + N)./max(Lpop,1e-12);

        %% OU mean paths
        beta1_t = p.beta1_bar + (p.beta10-p.beta1_bar)*exp(-p.theta_beta1*t);
        beta2_t = p.beta2_bar + (p.beta20-p.beta2_bar)*exp(-p.theta_beta2*t);
        beta3_t = p.beta3_bar + (p.beta30-p.beta3_bar)*exp(-p.theta_beta3*t);

        eta1_t = eta1_bar + (eta10-eta1_bar)*exp(-p.theta_eta1*t);
        eta2_t = p.eta2_bar + (p.eta20-p.eta2_bar)*exp(-p.theta_eta2*t);

        %% Brownian increments
        dWb1 = sqrt_dt*randn(1,nPaths);
        dWb2 = sqrt_dt*randn(1,nPaths);
        dWb3 = sqrt_dt*randn(1,nPaths);
        dWe1 = sqrt_dt*randn(1,nPaths);
        dWe2 = sqrt_dt*randn(1,nPaths);

        %% Drift
        drift = zeros(9,nPaths);
        diff  = zeros(9,nPaths);

        drift(1,:) = (1-p.rho)*p.Lambda_S.*(1-Ppop/p.K_S) - (p.alpha_2*A + p.gamma_1*lambda + beta2_t*A + p.kappa + p.mu_S).*SU;
        drift(2,:) = p.rho*p.Lambda_S.*(1-Ppop/p.K_S) + p.kappa*SU - (p.alpha_1*A + beta3_t*A + p.gamma_2*lambda + p.mu_S).*SI;
        drift(3,:) = (p.alpha_1*SI + p.alpha_2*SU).*A - (p.gamma_3*lambda + beta1_t*A + p.mu_S).*C;
        drift(4,:) = (beta1_t*C + beta2_t*SU + beta3_t*SI).*A - p.mu_S*R;
        drift(5,:) = (p.gamma_1*SU + p.gamma_2*SI + p.gamma_3*C).*lambda;
        drift(6,:) = p.Lambda_L.*(1-Lpop/p.K_L) - (p.pi+p.mu_L).*T;
        drift(7,:) = (1-p.omega)*p.pi*T + p.nu_2*U + p.nu_1*N - (p.epsilon+p.mu_L+eta1_t*lambda).*A;
        drift(8,:) = p.omega*p.pi*T + p.epsilon*A  - (p.nu_2+p.mu_L+eta2_t).*U;
        drift(9,:) = eta1_t*lambda.*A + eta2_t*U  - (p.nu_1+p.mu_L).*N;

        %% Diffusion terms

        % beta_2 noise: SU -> R
        diff(1,:) = diff(1,:) - phi_beta2(k).*A.*SU.*dWb2;
        diff(4,:) = diff(4,:) + phi_beta2(k).*A.*SU.*dWb2;

        % beta_3 noise: SI -> R
        diff(2,:) = diff(2,:) - phi_beta3(k).*A.*SI.*dWb3;
        diff(4,:) = diff(4,:) + phi_beta3(k).*A.*SI.*dWb3;

        % beta_1 noise: C -> R
        diff(3,:) = diff(3,:) - phi_beta1(k).*A.*C.*dWb1;
        diff(4,:) = diff(4,:) + phi_beta1(k).*A.*C.*dWb1;

        % eta_1 noise: A -> N
        diff(7,:) = diff(7,:) - phi_eta1(k).*lambda.*A.*dWe1;
        diff(9,:) = diff(9,:) + phi_eta1(k).*lambda.*A.*dWe1;

        % eta_2 noise: U -> N
        diff(8,:) = diff(8,:) - phi_eta2(k).*U.*dWe2;
        diff(9,:) = diff(9,:) + phi_eta2(k).*U.*dWe2;

        %% Euler-Maruyama update
        Y = Y + drift*dt + diff;
        Y = max(Y,0);

    end

    Nfinal = Y(9,:);
    Ufinal = Y(8,:);
    Dfinal = Ufinal + Nfinal;
    Psurv(j) = mean(Dfinal > Ntol);
    Nmean(j) = mean(Nfinal);
    Dmean(j) = mean(Dfinal);

    fprintf('eta1 %.3f | R0s %.3f | Psurv %.3f | mean N %.4f\n', eta1_bar, R0s_vals(j), Psurv(j), Nmean(j));
end

%% Sort for plotting
[R0s_sort, idx] = sort(R0s_vals);
Psurv_sort = Psurv(idx);
%Nmean_sort = Nmean(idx);
%eta1_sort = eta1_vals(idx);
Dmean_sort = Dmean(idx);

%% survival probability vs eta_1
figure;
plot(eta1_vals, Psurv, 'r-', 'LineWidth', 5)
xlabel('\eta_1')
ylabel('Discouragement survival probability')
grid on; box on;
set(gca,'FontSize',12,'LineWidth',2)

%% Survival probability vs discouragement population
figure;
scatter(Dmean_sort,Psurv_sort, 100,'filled')
hold on
plot(Dmean_sort, smoothdata(Psurv_sort,'movmean',5), 'LineWidth',4)
xlabel('Discouragement population, (U(t) + N(t))')
ylabel('Discouragement survival probability')
legend('Simulated data', 'Probability', Location='best')
grid on
box on
xlim([0 100])
set(gca,'LineWidth',2)

end