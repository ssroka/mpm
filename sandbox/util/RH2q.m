
function  [q] = RH2q(RH,T,P)


% reference : Pruppacher and Klett (1978) [abbrv. PK78]
% 
% sanity check : http://go.vaisala.com/humiditycalculator/5.0/
%
% input: 
%   RH = relative humidity as a percent (0<= RH <= 100)
%    T = temperature in Celsius 
%    P = pressure in [mb]
% 
% output:
%   q = specific humidity [non-dim]
%
%

% PK78 eq (4-34)
% RH    = w_v / w_v_sat = e / e_sat
% w_v   = mixing ratio
% w_v   = saturation mixing ratio
% e     = water vapor pressure
% e_sat = saturation water vapor pressure


% PK78 eq (4-24)
% w_v = rho_v / rho_a
% rho_v = water vapor density
% rho_a = dry air density

% Ideal Gas Law : e = rho_v R_v T


%% Step 1: calculate e_sat for water and ice in [mb]
% Buck 1981 equation for e_sat

e_sat_H2O = (1.0007+3.46e-6*P)*(6.1121*exp(17.502*T/(240.97+T)));
e_sat_ICE = (1.0003+4.18e-6*P)*(6.1115*exp(22.452*T/(272.55+T)));

% by convention, not in Buck 1981
if T>=0
	e_sat_tot = e_sat_H2O;
elseif T<=i-30
	e_sat_tot = e_sat_ICE;
else % linear combination
	e_sat_tot = e_sat_H2O*(T/30 + 1) + e_sat_ICE*(-T/30)
end

%% Step 2: calculate e in [mb]

% PK78 eq (4-34)
e = RH/100*e_sat_tot;


%% Step 3: calculate w_v in [none]

M_a = 0.0289644; % [kg/mol] molecular weight of water; Andreas 2005 Symbols
M_w = 0.018015;  % [kg/mol] molecular weight of water; Andreas 2005 Symbols
% PK78 near eq (4-30)
epsilon = M_w/M_a;

% PK78 eq (4-40)
w_v = epsilon*e/(e+P)


%% Step 4: calculate q in [none]

% PK78 eq (4-38)
q = w_v/(1+w_v);







end









































