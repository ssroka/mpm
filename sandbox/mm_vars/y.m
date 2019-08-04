function [y] = y(T_s,T_a,r,m_s,s,P)

% function [y] = y(T_s,T_a,r,m_s,s,P)
%  y = models how curvature (the Kelvin effect) and dissolved salt affect
%      the saturation vapor pressure at the surface of an aqueous solution 
%      droplet
%
% reference : (Andreas 2005, eq 2.2)
%
% input: 
%     T_s = sea water temperature [deg C]
%     T_a = ambient air temperature [deg C]
% sigma_s = surface tension [N/m]                         
%           ** function call to sigma_s.m, input T_s [deg C], m_s [kg],
%                                           m_w [kg]
%  rho_w  = density of pure water [kg/m^3]                 
%             ** function call to rho_w.m, input is T [deg C] and P [Pa].
%      r  = instantaneous droplet radius in [m]              
%   Phi_s = the practical osmotic coefficient of NaCl 
%           dissolved in water
%           ** function call to Phi_s.m, input is m [moles/kg]
%     m_s = mass of the dissolved salt (constant) [kg]
%       s = S/1000 for S salinity in ppt
%   rho_s = density of the saline drop [kg/m^3]            
%             ** function call to rho_s.m, input is T [deg C], r [m],
%                                                   m_s [kg], s [kg/kg]
%       P = pressure [Pa]
%
% constants:
%     M_w = 0.018015 molecular weight of water [kg/mol]   
%      R  = 8.31441 universal gas constant [J/(mol K)]    
%     nu  = 2, the number of ions NaCl dissociates into   
%     M_s = 0.058443 molecular weight of salt [kg/mol]   
%             * loaded from microphysicalConstants.mat
%
% output:
%       y = [dimensionless]
%

load('microphysicalConstants.mat','M_w','R','nu','M_s')

T_a_Kelvin = T_a+273.15; %[K]

m_w = m_s*(1-s)/s;

% molality is needed for Phi_s. NOTE there is a typo in Andreas 2005 eq 2.4
m = m_s/(M_s*m_w);

y = 2.*M_w.*sigma_s(T_s,m_s,m_w)./(R.*T_a_Kelvin.*rho_w(T_s,P).*r)-(nu.*Phi_s(m).*m_s.*M_w./M_s)./((4.*pi.*rho_s(T_s,r,m_s,s,P).*r.^3./3)-m_s);

end
