function [Phi_s] = Phi_s(m)
 % function [Phi_s] = Phi_s(m)
% See Andreas 1989 equation 27
%
%
% input:
%      m  = molality [moles of salt / kg of water]
%
% output:
%   Phi_s = practicle osmotic coefficient  [kg / m^3]
%
%

Phi_s = 0.927-2.164*(10^(-2))*m+3.486*(10^(-2))*m^2-5.956*(10^(-3))*m^3+3.911*(10^(-4))*m^4;
end
