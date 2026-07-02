function R=RRR(LHSmatrix,x)

Parameter_settings_LHS;
pi=LHSmatrix(x,1);
nu_1=LHSmatrix(x,2);
eta_1=LHSmatrix(x,3);
mu_L=LHSmatrix(x,4);
sigma_1=LHSmatrix(x,5);
theta_1=LHSmatrix(x,6);

Rs=(eta_1*pi)/((nu_1+mu_L)*(pi+mu_L)) - (sigma_1^2*pi^2)/(4*theta_1*(nu_1+mu_L)*(pi+mu_L)^2);

R=(Rs);
end