function [rho] = rho_s(T,r,m_s,s,P)
% function [rho] = rho_s(T,r,m_s,s,P)
% 
% input: 
%     T     = sea water density    [deg C]
%     r     = radius of the drop   [m]
%     m_s   = mass of salt in drop [kg]
%     s     = salinity of drop     [kg/kg]
%     rho_w = density of pure water [kg m^-3]
%             .*.* function call to rho_w.m, input is T [deg C]
%     P     = ambient pressure     [Pa]
%
% constants:
%     M_s = molecular mass of NaCl [kg./mol]
%             .* loaded from microphysicalConstants.mat
%
% output: 
%     rho = density of sea water [kg m^-3]
global Nayar_flag;
if Nayar_flag
	rho = SW_Density(T,'C',s,'w',P,'Pa');
else
	% reference Andreas 1989 , eq 12 to eq 17
	load('microphysicalConstants.mat','M_s')

	% -----------------------------------------
	% from Andreas 1989 -----------------------
	v_a0 = 1e-6.*(12.97+0.2340.*T-4.210.*10^(-3).*T.^2 + 2.857.*10^-5 .* T.^3); %eq 16
	Sv   = 1e-6.*( 2.982 - 4.970 .* 10^(-2).*T + 6.032 .* 10^(-4).*T.^2); %eq 17
	c    = 10^(-3) .* m_s./M_s./(4.*pi.*r.^3./3); % eq 15
	v_a  = v_a0 + Sv.*sqrt(c); % eq 14

	% -----------------------------------------
	% from Andreas 2005 -----------------------
	m_w = m_s.*(1-s)./s;  % eq 2.6

	rho = rho_w(T) .* (1+m_s./m_w)./(1+v_a.*rho_w(T)./M_s.*m_s./m_w); %eq 12

end










end
