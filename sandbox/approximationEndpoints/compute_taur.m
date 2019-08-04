% --- compute tau_r ---------------------------------------------------
%
% reference: Andreas 2005, eq 6.8
%
% input:
%           f = fractional RH  RH/100             **[Teq.m]
%           y = the Kelvin effect [none]                    **[Teq.m]
%         L_v = latent heat of vaporization [J kg^-1]       **[mC.m]
%           r = instantaneous droplet radius
%  k_a_prime = thermal conductivity of water accounting     **[Teq.m]
%              for non-continuum effects [W/(m K)]
%        drdt = the time-derivative of the droplet radius given in Andreas
%               2005, eq 2.1
%       d2rdt2 = the second derivative of the droplet radius evolution with
%                time, easily computed from drdt (note we use drdt to get d2rdt2)
%         r_0 = initial droplet radius [m]                 **[in.m]
%        r_eq = equilibrium droplet radius [m]                     
%
% output:
%  tau_r = the time between ejection and when the equilibrium radius is reached [s]
%
T_a_Kelvin = Celsius2Kelvin(T_a); %[K]
P_mb       = Pa2mb(p0);
e_sat_Pa   = mb2Pa(e_sat(T_a,P_mb));

load('microphysicalConstants.mat','M_w','R')

s =@(r)m_s/(m_s+compute_m_w(T_s_0,r,m_s,p0,maxEr_s,maxIt));

drdt_numerator   =@(r) (((f-1)-y(T_s_0,T_a,r,m_s,s(r),p0))./r);
drdt_denominator =@(r) (rho_s(T_s_0,r,m_s,s(r),p0).*R.*T_a_Kelvin./(D_w_prime(T_s_0,P_mb,r).*M_w.*e_sat_Pa)+rho_s(T_s_0,r,m_s,s(r),p0).*L_v(T_s_0,s(r)*1000)./(k_a_prime(T_s_0,P_mb,r).*T_a_Kelvin).*(L_v(T_eq-273,s(r)*1000).*M_w./(R.*T_a_Kelvin)-1));
drdt             =@(r) drdt_numerator(r)./drdt_denominator(r);

integrate_drdt;

ele_to_fit = min(130,length(r_vec));%max(round(length(time_vec_r)*0.1),10);
orderOfFit=min(4,length(r_vec)-1);
f2 = polyfit(time_vec_r(1:ele_to_fit),r_vec(1:ele_to_fit),orderOfFit);%polyfit(time_vec_r(1:end),r_vec(1:end),orderOfFit);
d2rdt2 = 2*f2(orderOfFit-1);
%{
semilogx(time_vec(1:ele_to_fit),sum([time_vec(1:ele_to_fit).^4 time_vec(1:ele_to_fit).^3 time_vec(1:ele_to_fit).^2 time_vec(1:ele_to_fit) ones(ele_to_fit,1)]*f2',2),'gs')
%}

tau_r = (-drdt(r_0)-(3*drdt(r_0).^2-2*(r_0-r_eq).^2*d2rdt2).^(0.5))/(d2rdt2-(drdt(r_0).^2/(r_0-r_eq)));

