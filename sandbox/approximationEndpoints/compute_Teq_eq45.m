T_eq_vec = zeros(size(T_s_t));
T_eq_vec(1) = T_a+L_v(T_s)*D_w_prime(T_s,P_mb,r_0)*(rho_v-rho_vr)/k_a_prime(T_s,P_mb,r_0);

e_sat_Pa = mb2Pa(e_sat(T_a,P_mb));
rho_v  = (M_w*e_sat_Pa/(R*T_a_Kelvin))*f;

for i = 1:length(T_s_t)-1
    rho_vr = (M_w*mb2Pa(e_sat(T_eq_vec(i),P_mb))/(R*Celsius2Kelvin(T_eq_vec(i))))*exp(y(T_eq_vec(i),T_a,r_t(i+1),m_s,s(r_0),ic.p0));
    T_eq_vec(i+1) = 0.5*(T_eq_vec(i) + T_a + L_v(T_s_t(i))*D_w_prime(T_s_t(i),P_mb,r_t(i))*(rho_v-rho_vr)/k_a_prime(T_s_t(i),P_mb,r_t(i)));
end

kk = find(time_vec<tau_T);
T_eq_vec(max(kk)+1)
KK = find(time_vec<tau_r);







