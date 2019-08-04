
function [L_v] = L_v(T,S)
% function [L_v] = L_v(T,S)
%
%
% input:
%      T = sea water temperature in degrees Celsius
%      S = salt content in [ppt]
% output:
%   L_v = latent heat of vaporization in J/kg
%
global Nayar_flag;

if Nayar_flag
    L_v = SW_LatentHeat(T,'C',S,'ppt');
else
    % reference : (Andreas 1989, equation 35)
    L_v = (25.00-0.02274.*T).*10^5; % [J /kg]
end
end
