function [k_a_prime] = k_a_prime(T,P,r)

% function [k_a_prime] = k_a_prime(T,P,r)
%  k_a_prime = thermal conductivity of air accounting for non-continuum
%  effects [W/(m K)]
%
% reference : (Andreas 1995, eq 4.3)
%
% empirical constants:
%     Delta_T = 2.16*(1e-7); % [m]
%     alpha_T = 0.7;         % [none]
%
% input: 
%    k_a = thermal conductivity [W/(m K)] of water as a function of ambient 
%          temperature [deg Celsius] Andrea 1995 eq 4.5
%     r  = instantaneous droplet radius in m              
%  rho_a = density of dry air [kg/m^3]                   
%             ** function call to rho_a.m, input is T [K] and P [Pa]
%      T = sea water temperature in degrees Celsius  
%      P = pressure in [mb]
%
% constants:
%     R  = 8.31441 universal gas constant in J/(mol K)    
%    M_a = 0.0289644 molecular weight of air in [kg/mol] 
%             * loaded from microphysicalConstants.mat
%   c_pd = 1006 specific heat of air in J/(kg K)
% 
% output:
%   k_a_prime = saturated vapor pressure in W/(m K)
%
helpingAnonFxns;

Delta_T = 2.16*(1e-7); % [m]
alpha_T = 0.7;         % [none]

T_s_Kelvin = T + 273.15; %[K]

load('microphysicalConstants.mat','M_a','R','c_pd')

P_Pa = mb2Pa(P); % convert mb to Pa

k_a       = @(T) 0.02411*(1 + 3.309*1e-3*T-1.441*1e-6*T^2); % [W/(m K)]
k_a_prime = k_a(T)./(r./(r+Delta_T)+(k_a(T)./(r.*alpha_T.*rho_a(T_s_Kelvin,P_Pa).*c_pd)).*sqrt(2.*pi.*M_a./(R.*T_s_Kelvin)));



end