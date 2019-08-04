function [tau_f] = compute_tauf(U10,T_s_0,r_0,m_s,s0,p0,T_a,maxEr_uf,maxIt)
% --- compute tauf ---------------------------------------------------
% [tau_f] = compute_tauf(U10,T_s_0,r_0,m_s,s0,p0,T_a,maxEr_uf,maxIt)
%
% reference: Andreas 1992, eq 5
%            
% input:
%    A13 = Characteristic wave amplitude [m]          
%    U10 = 10-m windspeed [m/s]                       
%    u_f = Stokes speed of fall for larger Reynolds Numbers [m/s] 
%          Andreas 1990, eq 5.1: u_f = (2*r_0^2*g)/(9*nu_a*(1+0.158*(2*r_0*u_f/nu_a)^(2/3)))*(rho_s/rho_a - 1);
%          Andreas 1989, eq 81
%    r_0 = initial droplet radius [m]      
%  rho_s = density of seawater [kg m^-3]
%     p0 = ambient pressure in [Pa]
%   nu_a = the kinematic viscosity of air [m^2/s] evaulated at air 
%          temperature according to Andreas 1989 eq 78  
%    s0  = salinity in [kg/kg]
%
% constants:
%      g = gravity 9.81 [m/s^2]
%
% output:
%  tau_f = the time of flight of the droplet in seconds
%
load('microphysicalConstants.mat','g')
helpingAnonFxns;

A13 = 0.015*U10^2; % [m]

% might need density of moist air NOT density of dry air
T_s_Kelvin = Celsius2Kelvin(T_s_0);

fofuf =@(u_f) u_f - (2*g*r_0^2*(rho_s(T_s_0,r_0,m_s,s0,p0)/rho_a(T_s_Kelvin,p0) - 1))/(9*nu_a(T_a)*((0.158*((2*r_0*u_f)/nu_a(T_a))^(2/3)) + 1));
dfduf =@(u_f) (158*g*r_0^3*(rho_s(T_s_0,r_0,m_s,s0,p0)/rho_a(T_s_Kelvin,p0) - 1))/(3375*nu_a(T_a)^2*((0.158*((2*r_0*u_f)/nu_a(T_a))^(2/3)) + 1)^2*((2*r_0*u_f)/nu_a(T_a))^(1/3)) + 1;

% Newton's method
ItNum = 1;
u_f_old = calcVterm_sphere(r_0,rho_s(T_s_0,r_0,m_s,s0,p0),g);
error_uf = abs(fofuf(u_f_old));
while error_uf > maxEr_uf/1000 && ItNum < maxIt
    u_f = u_f_old - fofuf(u_f_old)/dfduf(u_f_old);
    error_uf = abs(fofuf(u_f));
    u_f_old = u_f;
    ItNum = ItNum + 1;
end
u_f = u_f_old;


tau_f = A13/u_f;









