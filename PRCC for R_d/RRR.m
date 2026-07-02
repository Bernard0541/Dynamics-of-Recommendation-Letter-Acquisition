function R=RRR(LHSmatrix,x)

Parameter_settings_LHS;
pi=LHSmatrix(x,1);
nu_1=LHSmatrix(x,2);
eta_1=LHSmatrix(x,3);
mu_L=LHSmatrix(x,4);

Rd=(eta_1*pi)/((nu_1+mu_L)*(pi+mu_L));

R=(Rd);
end