function specificHeatMoistAir = c_pm(T_a_vec,P, RH)
% function specificHeatMoistAir = c_pm(T_a_vec,P, RH)
%
% input:
%     T_a   = moist air temperature    [deg C]
%     P     = total pressure in [Pa]
%     RH    = relative humidity [%]
%c_pm(
% constants:
%     R     = universal gas constant
%
% output:
%     specificHeatMoistAir = specific heat of moist air  [J kg^-1 K^-1]
%
% reference Atmospheric Convection, ch4 Moist Thermodynamics and Stability
%           Kerry Emanuel

%

load('microphysicalConstants.mat','c_pd');

if logical(sum(double(T_a_vec<0))) || logical(sum(double(T_a_vec>30)))
    warning('\nSaturated vapor pressure calculation\n is valid for temperatures between 0 and 30 degrees C\n')
end

R_v = 461.51; % J./ (kg K) from Atmospheric Convection
run ./util/helpingAnonFxns.m
P_mb = Pa2mb(P);
specificHeatMoistAir = zeros(length(T_a_vec),1);
for jj = 1:length(T_a_vec)
    T_a = T_a_vec(jj);
    T_a_Kelvin = Celsius2Kelvin(T_a);
    e = (RH(jj)./100).*e_sat(T_a,P_mb(jj));
    rho_v = mb2Pa(e)./(R_v.*T_a_Kelvin); % vapor density
    rho_d = rho_a(T_a_Kelvin,P(jj)); % dry air density
    r = rho_v./rho_d;
    specificHeatMoistAir(jj) = c_pd.*(1+0.85.*r);
end

end
































