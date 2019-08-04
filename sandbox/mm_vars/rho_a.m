function [rho_a] = rho_a(T,P)
% function [rho_a] = rho_a(T,P) density of dry air
% See Andreas 1989 equation 31
%
%   The density of dry air is computed for different temperatures
%   and pressures.
%
% input:
%    T    = temperature in Kelvin
%    P    = pressure in Pa
%
% constants:
%    T0   = reference temperature in Kelvin
%    P0   = reference pressure in Pa
%
% output:
%    rho_a = density  [kg / m^3]
%
%
P0 = 101325; % Pa
T0 = 273.15; % K
rho_a = 1.2923*(T0/T)*(P/P0);
