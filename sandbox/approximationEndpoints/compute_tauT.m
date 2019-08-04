% --- compute tauT ---------------------------------------------------
%
% reference: Andreas 2005, eq 4.11
%
% empirical constants:
%           a = 17.502; % [none] 
%           b = 240.97; % [deg Celsius]
%
% input:
%       rho_s = droplet density [kg/m^3]                          
%         r_0 = initial droplet radius [m]                
%   k_a_prime = thermal conductivity of water accounting for non-continuum
%                effects [W/(m K)]                         
%         L_v = latent heat of vaporization               
%   D_w_prime = diffusivity of water accounting for non-continuum
%               effects [m^2/s]                           
% drho_vsatdt = the time rate of change in saturation vapor density [kg/(m^3 s)]
%                Andreas 2005 eq 4.8
%    rho_vsat = the saturation vapor density [kg/m^3] Andreas 2005 eq 2.11
%           y = models how curvature (the Kelvin effect) and dissolved salt affect
%               the saturation vapor pressure at the surface of an aqueous solution 
%               droplet         
%
% constants:
%      R  = 8.31441 universal gas constant in J/(mol K)
%     M_w = molecular weight of water [kg/mol]                     
%   
%
% output:
%  tau_T = the time between ejection and when the droplet reaches thermal
%  equilibrium [s]
%


load('microphysicalConstants.mat','M_w','R')

T_a_Kelvin = Celsius2Kelvin(T_a);
T_s_Kelvin = Celsius2Kelvin(T_s);
P_mb       = Pa2mb(p0);

e_sat_Pa = mb2Pa(e_sat(T_a,P_mb));
rho_vsat = (M_w*e_sat_Pa/(R*T_a_Kelvin))*exp(0);% where y = 0

drho_vsatdT = rho_vsat*((a*b/((b+T_a_Kelvin-273.15)^2)) - 1/(T_a_Kelvin)); % T_a in K

tau_T = (rho_s(T_s,r_0,m_s,s0,p0)*c_ps(T_s,S0,p0)*r_0^2)/(3*(k_a_prime(T_s,P_mb,r_0)+L_v(T_s,S0)*D_w_prime(T_s,P_mb,r_0)*drho_vsatdT));











