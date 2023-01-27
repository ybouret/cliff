-- pH time scale
A_h = 24.7; -- s
B_h = 0.034; -- /mM

-- pH amplitude
pH_ini   = 5.92;
pH_ref   = 7.40; -- set to pH_ini to fix pH
Lambda_h = 15.2; -- mM

-- solution
rho_s        = 12.0192;
-- sigma        = 1.0/0.99772;
sigma        = 1.003772;
sigma        = 1.00377242943013
d7out        = 14.57;
Temperature  = 273.5+37;

-- integration
ftol = 1e-6;
ctrl = 0.1;

-- outer pH
pH_out = 7.4
