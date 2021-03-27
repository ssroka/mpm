
saveDir = 'Exp7220';

% print the error between the calculated microphysical quantity and the
% quantity stated in the corresponding publication
print_flag = true;

Nayar_flag = true;

% SSGF String
SSGF_str = 'singleDrop';

RH = 82.0;

% ambient air temperature
T_a = 27.5;

% initial droplet temperature (and local SST)
T_s_0 = 29.0;

% initial salinity of drop
S0 = 34; % ppt 

% pressure of air
p0 = 100000; % Pa

% initial droplet radius
r_0_vec = [50:-50:50]*1e-6; % m

% 10-m wind speed
U10 = 20.0; % m/s

% integration time
t_final = 2050; % s

% Newton-Raphson for r_eq
maxIt = 1000; % the maxium number of iterations before aborting the iterative scheme

maxEr_req = max(r_0_vec)*0.00001; % the max error in r_eq
maxEr_uf  = calcVterm_sphere(max(r_0_vec))*0.0001;  % the max error in u_f
maxEr_s   = S0/1000*0.01; % the max error

 % micro-physical endpoints from caption in Fig 11 of Andreas 2005
T_eq_exact = 17.07; % deg C
r_eq_exact = 61.44* 10^-6; % m
tau_r_exact = 303.0; % s
tau_T_exact = 0.0176; % s