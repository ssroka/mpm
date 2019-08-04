function [esat] = e_sat(T,P)
% function [esat] = e_sat(T,P)
% reference : Andreas 1995 eq 2.3
%
% empirical constants:
%     a = 17.502; % [none] 
%     b = 240.97; % [deg Celsius]
%
% input: 
%    T =  sea or air temperature in Celsius 
%    P = pressure in [mb]
%
% output:
%   esat = saturated vapor pressure in [mb]
%
a = 17.502; % [none] 
b = 240.97; % [deg Celsius] 

% Buck 1981, eq 2a-b, 8, where T_a is in deg C and P is in mb
esat = (1.0007 + 3.46*(1e-6)*P) * 6.1121*exp((a*T)/(b+T));
% for some reason Andreas 1995 eq 2.3 sets e_sat = 6.1375exp((a*T_a)/(b+T_a))
%e_sat =@(T_a) 6.1375*exp((a*T_a)/(b+T_a));
