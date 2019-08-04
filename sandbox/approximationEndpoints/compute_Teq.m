% --- compute T_eq ---------------------------------------------------
%
% input:
%    T_a = ambient air temperature in degrees Celsius
%     p0 = ambient air pressure in Pa
%     r  = instantaneous droplet radius in m
%
% output:
%   T_eq = droplet equilibrium temperature in Kelvin
%

T_a_Kelvin = Celsius2Kelvin(T_a); %[K]
T_s_Kelvin = Celsius2Kelvin(T_s); %[K]
P_mb       = Pa2mb(p0);           %[mb]

%% compute beta

% --- compute beta ---------------------------------------------------
% reference: Andreas 2005, eq 3.3
e_sat_Pa = mb2Pa(e_sat(T_a,P_mb));

load('microphysicalConstants.mat','M_w','R')
beta  = e_sat_Pa*L_v(T_s,S0)*M_w*D_w_prime(T_s,P_mb,r_0)/(T_a_Kelvin*R*k_a_prime(T_s,P_mb,r_0));
clear M_w R

%
% reference: Andreas 2005, eq 3.1
%
% empirical constants:
%     a = 17.502; % [none] 
%     b = 240.97; % [deg Celsius]
%
% input:
%   beta = [K]
%    T_a = ambient air temperature in K
%  alpha = [none]
%      y = the Kelvin effect [none]
%      f = fractional RH (see in.m) RH/100
%
% output:
%   T_eq = droplet equilibrium temperature in Kelvin
%

f = RH/100; %[none]

a = 17.502; % [none] 
b = 240.97; % [deg Celsius] 
alpha = a*b*T_a_Kelvin/((b+T_a_Kelvin-273.15)^2);

T_eq_a = (beta/(T_a_Kelvin^2))*(alpha^2/2-alpha*((2*T_a_Kelvin+b-273.15)/(T_a_Kelvin+b-273.15))+1)*exp(y(T_s,T_a,r_0,m_s,s0,p0));
T_eq_b = 1+beta/T_a_Kelvin*(alpha-1)*exp(y(T_s,T_a,r_0,m_s,s0,p0));
T_eq_c = -beta*(f-exp(y(T_s,T_a,r_0,m_s,s0,p0)));


quadFormula =@(a,b,c) (-[b;b] + [sqrt(b^2-4*a*c);-sqrt(b^2-4*a*c)])./(2.*[a;a]);

DeltaT = quadFormula(T_eq_a,T_eq_b,T_eq_c);

T_eq_tmp = DeltaT+T_a_Kelvin;
T_eq     = T_eq_tmp((T_eq_tmp-273.15)>0);% remove the negative temperature ('reminder we are in Kelvin')

if sum(double(((T_eq_tmp-273.15)<0)))==2 || isempty(T_eq)
    error('No positive T_eq found')
end





