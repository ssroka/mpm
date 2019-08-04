function [D_w_prime] = D_w_prime(T,P,r)
%  function [diffusivity] = D_w_prime()
%  D_w_prime = diffusivity of water accounting for non-continuum
%  effects [m^2/s]
% 
% reference : (Andreas 1995, eq 4.2)
%             (Andreas 2005, pg 332 under eq 4.7 indicates
%              that the arguments are r_0 and T_s)
%
% empirical constants:
%     Delta_w = 8*(1e-8); % [m]
%     alpha_c = 0.036;    % [none]
%         T0  = 273.15, reference temperature [K]
%         P0  = 1013.25, reference pressure [mb or hPa]
%
% input: 
%    P = ambient pressure [mb] 
%    r = instantaneous droplet radius   [m]              
%    T = sea water temperature in [degrees Celsius]
%
% constants:
%   R  = 8.31441 universal gas constant in J/(mol K)    
%  M_w = 0.018016 molecular weight of water in [kg/mol] 
%             * loaded from microphysicalConstants.mat
%
% output:
%   D_w_prime = saturated vapor pressure in m^2/s
%

Delta_w = 8*(1e-8); % [m]
alpha_c = 0.036;    % [none]
T0  = 273.15; % [K]
P0  = 1013.25;% [mb]

T_s_Kelvin = T + 273.15; %[K]

load('microphysicalConstants.mat','M_w','R')

D_w         =@(T,P) 2.11*(1e-5)*((T+273.15)/T0)^1.94*P0/P;% Andreas 1989  eq 29
D_w_prime   =D_w(T,P)./(r./(r+Delta_w)+(D_w(T,P)./(r.*alpha_c)).*sqrt(2.*pi.*M_w./(R.*T_s_Kelvin)));



end