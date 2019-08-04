function [sigma_s] = sigma_s(T_s,m_s,m_w,s)
% reference: Andreas 1989
%
% input:
%   T_s = temperature of drop in deg Celsius
%   m_s = mass of salt in drop in kg
%   m_w = mass of water in drop in kg
%  (optional)
%   s   = salinity of water in kg/kg = m_s/(m_s+m_w)
%
% output:
%   sigma_s = surface tension in N/m
%
% NOTE!! - I believe there is a typo in Andreas 1989 eq 33
% I looked at Pruppacher and Klett 1978, p. 107 and when I convert between
% cgs and SI I get 2.77*1e-2 NOT 2.77*1e-5


% if T_s>40 || T_s<-40
%     error('surface tension not coded for temperatures outside -40 deg C\n to 40 deg C')
% end
global Nayar_flag

if Nayar_flag
    if nargin < 4
        sigma_s = SW_SurfaceTension(T_s,'C',m_s/(m_s+m_w),'w')/1000; % [mN/m]/1000 = [N/m]
    else
        sigma_s = SW_SurfaceTension(T_s,'C',s,'w')/1000; % [mN/m]/1000 = [N/m]
    end
else
    sigma_w = 7.610*(1e-2)-1.55*(1e-4)*T_s;
    sigma_s = sigma_w + 2.77*1e-2*(m_s/m_w);
end








