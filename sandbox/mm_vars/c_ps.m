function specificHeat = c_ps(T_s,S,P)
% function specificHeat = c_ps(T_s,S,P)
% 
% input: 
%     T   = water temperature    [deg C]
%     S   = salt content of sea water in [ppt]
%     P   = pressure in [Pa]
%
% output: 
%     specific heat of sea water  [J kg^-1 K^-1]

global Nayar_flag

if Nayar_flag
	specificHeat = SW_SpcHeat(T_s,'C', S,'ppt',P,'Pa'); % [J kg^-1 K^-1]
else
	% page v of Andreas 1989
	specificHeat = 4*10^3; % [J kg^-1 deg C ^-1]

end







end
