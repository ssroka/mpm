function [rho] = rho_w(T,P)
% function [rho] = rho_w(T,P)
% 
% input: 
%     T   = water temperature    [deg C]
%     P   = pressure in [Pa]
%
% output: 
%     rho = density of pure water [kg m^-3]

global Nayar_flag

if Nayar_flag
	rho = SW_Density(T,'C',0,'ppt',P,'Pa');
else
	% reference Andreas 1989 , eq 13a
	rho = (999.8396 + 18.224944 * T- 7.922210 * 10^(-3)* T.^2)./(1+1.8159725*10^(-2).*T); %eq 13a
end











end
