function [nu_a] = nu_a(T_a)
% 
% reference : Andreas 1989 eq 78
%
%
% input: 
%    T =  air temperature in Celsius 
%
% output:
%   nu_a = viscosity of air in m^2 /s
%

if T_a < -173 || T_a > 277
    error('\n\n viscosity not valid for temperature T_a = %f',T_a)
end


nu_a = 1.326*10^-5*(1+6.542*1e-3*T_a+8.301*1e-6*T_a^2-4.840*1e-9*T_a^3);
