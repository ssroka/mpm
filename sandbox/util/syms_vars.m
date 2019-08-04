% all syms
clear; clc
syms M_w sigma_s R T_a_Kelvin T_s_Kelvin rho_wp r_eq nu Phi_s m_s M_a M_s rho_s r_eq rho_a PI c_pd S
syms f y r rho_s  T_a D_w_prime M_w e_sat e_sat_Pa L_v k_a_prime t k_a D_w Delta_T Delta_w alpha_T alpha_c
syms u_f r_0 g nu_a rho_s rho_a a b c_ps tau_T T_s p0 U10 tau_r
syms drho_vsatdt rho_vsat dFdr0 T_s_t r_t T_eq tau_f beta0 beta1 d2rdt2_syms drdt
beta = [beta0 ; beta1];
% TODO: make rho_s a function of time... I'm not sure how to do this right now

%% -------------------------- k_a_prime --------------------------------------
k_a_prime_2 = k_a/(r/(r+Delta_T)+(k_a/(r*alpha_T*rho_a*c_pd))*sqrt(2*PI*M_a/(R*T_s_Kelvin)));

%% -------------------------- D_w_prime_2 --------------------------------------
D_w_prime_2 = D_w/(r/(r+Delta_w)+(D_w/(r*alpha_c))*sqrt(2*PI*M_w/(R*T_s_Kelvin)));

%% -------------------------- y --------------------------------------
y_2 = 2.*M_w.*sigma_s./(R.*T_a_Kelvin.*rho_wp.*r)-(nu.*Phi_s.*m_s.*M_w./M_s)./((4.*pi.*rho_s.*r.^3./3)-m_s);

%% -------------------------- drdt --------------------------------------
%drdt =  (((f-1)-y)./r)/((rho_s.*R.*T_a_Kelvin./(D_w_prime.*M_w.*e_sat_Pa)+rho_s.*L_v./(k_a_prime.*T_a_Kelvin).*(L_v.*M_w./(R.*T_a_Kelvin)-1)))

d2rdt2 = diff(drdt,r);



drdt_y = diff(drdt,y);
dydr = diff(y_2,r);

drdt_D_w_prime_2 = diff(drdt,D_w_prime);
d_D_w_prime_2_dr = diff(D_w_prime_2,r);

drdt_k_a_prime_2 = diff(drdt,k_a_prime);
d_k_a_prime_2_dr = diff(k_a_prime_2,r);

drdt_dr = diff(drdt,r);

d2rdt2 =(drdt_dr+drdt_y*dydr+drdt_D_w_prime_2*d_D_w_prime_2_dr+drdt_k_a_prime_2*d_k_a_prime_2_dr)*drdt 

tau_r = (-drdt-(3*drdt.^2-2*(r_0-r_eq).^2*d2rdt2_syms).^(0.5))/(d2rdt2_syms-1/(r_0-r_eq)*drdt.^2);





