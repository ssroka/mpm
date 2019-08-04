
% input file for Approximation forumlas for the microphysical properties of
% saline droplets
count = 1;
RH = 95; % relative humidity
    
% ambient air temperature
T_a = -10; % Celsius

% initial droplet temperature (and local SST)
T_s = 0; % Celsius

% salinity of drop 
S0 = 34; % psu

% pressure of air 
p0 = 100000; % Pa

% initial droplet radius
r_0_vec = 100*1e-6; % m

% 10-m wind speed
U10 = 20; % m/s


% time
t = 0.1; %s

% Newton-Raphson 
% for both u_f and r_eq
maxIt = 100; % the maxium number of iterations before aborting the iterative scheme

maxEr_req = max(r_0_vec)*0.0001; % the max error in r_eq
maxEr_uf  = calcVterm_sphere(max(r_0_vec))*0.001;  % the max error in u_f

% micro-physical endpoints from caption in Fig 11 of Andreas 2005
%T_eq_exact = 17.07;% deg C
%r_eq_exact = 61.44 * 10^-6;% m
tau_T_exact = 0.5;% s
tau_r_exact = 500;% s

