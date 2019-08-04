% density from Andreas 1989

function [m_s] = compute_initial_drop_conditions(T_s_0,S0,r_0,P)
% compute the initial density
%{
let   m_s/m_w0 = S/(1-S)
v_a = apparent molal volume
%}

global Nayar_flag

if Nayar_flag
    rho_s_0= SW_Density(T_s_0,'C',S0,'ppt',P,'Pa');
    m_s = 4/3*pi*r_0^3*rho_s_0*S0/1000;
else

s0 = S0/1000;
if T_s_0>40
    warning('This relation is only valid for T<= 40 degrees C,\n check Andreas 1989 for higher temperatures')
end
% initial mass of salt
load('microphysicalConstants.mat','M_s')
rho_w_0 = (999.8396 + 18.224944 * T_s_0- 7.922210 * 10^(-3)* T_s_0.^2)./(1+1.8159725*10^(-2).*T_s_0); %eq 13a
v_a0= 1e-6*(12.97+0.2340*T_s_0-4.210*10^(-3)*T_s_0.^2 + 2.857*10^-5 * T_s_0.^3); %eq 16
Sv=10^-6*( 2.982 - 4.970 * 10^(-2)*T_s_0 + 6.032 * 10^(-4)*T_s_0.^2); %eq 17
m_s = 4/3*pi*r_0^3*rho_w_0*s0/(1-s0); % eq 21
c = 10^(-3) * m_s/M_s/(4*pi*r_0^3/3); % eq 15
v_a = v_a0 + Sv*sqrt(c); % eq 14
rho_s0 = rho_w_0 * (1+s0/(1-s0))./(1+v_a*(rho_w_0./M_s)*s0./(1-s0)); %eq20 sub into eq 12
m_d0 = 4/3*pi*r_0^3*rho_s0;
m_s = s0*m_d0;
m_w_0 = m_d0 - m_s;
m_w = compute_m_w(T_s_0,r_0,m_s);
s = m_s/(m_w + m_s);
S = 1000*s;

end










