
% input file for Approximation forumlas for the microphysical properties of
% saline droplets

RH = 90; % relative humidity
    
% ambient air temperature
T_a = 25; % Celsius

% initial droplet temperature (and local SST)
T_s = 27; % Celsius

% salinity of drop 
S = 34; % psu

% pressure of air 
p0 = 101325; % Pa

% initial droplet radius
r_0 = 500*1e-6; % m

% 10-m wind speed
U10 = 20; % m/s

% SSGF
dFdr0 = 1/r_0; % 1/ (m^2 s^1 ) units

% time
t = 0.1; %s

% Newton-Raphson 
% for both u_f and r_eq
maxIt = 50; % the maxium number of iterations before aborting the iterative scheme

maxEr_req = r_0*0.001; % the max error in r_eq
maxEr_uf  = calcVterm_sphere(r_0)*0.001;  % the max error in u_f

% micro-physical endpoints from caption in Fig 11 of Andreas 2005
T_eq_exact = 17.07;% deg C
r_eq_exact = 61.44 * 10^-6;% m
tau_T_exact = 0.176;% s
tau_r_exact = 303;% s

