
% input file for Approximation forumlas for the microphysical properties of
% saline droplets

RH = 90; % relative humidity
    
% ambient air temperature
T_a = 18; % Celsius

% initial droplet temperature (and local SST)
T_s_0 = 20; % Celsius

% initial salinity of drop 
S0 = 34; % psu

% pressure of air 
p0 = 100000; % Pa

% initial droplet radius
r_0_vec = 100*1e-6;% m 

% integration time
t_final = 1050; % s

% Newton-Raphson for r_eq
maxIt = 1000; % the maxium number of iterations before aborting the iterative scheme

maxEr_req = max(r_0_vec)*0.00001; % the max error in r_eq

% micro-physical endpoints from caption in Fig 11 of Andreas 2005
T_eq_exact = 17.07;% deg C
r_eq_exact = 61.44 * 10^-6;% m
tau_T_exact = 0.176;% s
tau_r_exact = 303;% s

