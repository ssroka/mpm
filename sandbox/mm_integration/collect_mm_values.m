%% T_eq
min_T = min(stateVec(:,2));

T_eq_vec = zeros(size(T_s_t));
rho_v  = (M_w*e_sat_Pa/(R*T_a_Kelvin))*f;
rho_vr = (M_w*mb2Pa(e_sat(T_s,P_mb))/(R*Celsius2Kelvin(T_s)))*exp(y(T_s,T_a,r_0,m_s,s(r_0),ic.p0));

T_eq_vec(1) = T_a+L_v(T_s,s0)*D_w_prime(T_s,P_mb,r_0)*(rho_v-rho_vr)/k_a_prime(T_s,P_mb,r_0);
e_sat_Pa = mb2Pa(e_sat(T_a,P_mb));
rho_v  = (M_w*e_sat_Pa/(R*T_a_Kelvin))*f;
for i = 1:length(T_s_t)-1
    rho_vr = (M_w*mb2Pa(e_sat(T_eq_vec(i),P_mb))/(R*Celsius2Kelvin(T_eq_vec(i))))*exp(y(T_eq_vec(i),T_a,r_t(i+1),m_s,s(r_t(i+1)),ic.p0));
    T_eq_vec(i+1) = 0.5*(T_eq_vec(i) + T_a + L_v(T_s_t(i),s(r_t(i)))*D_w_prime(T_s_t(i),P_mb,r_t(i))*(rho_v-rho_vr)/k_a_prime(T_s_t(i),P_mb,r_t(i)));
end
kk = find(time_vec<tau_T);
T_eq_fmm = T_eq_vec(max(kk)+1);

%% r_eq
min_r = min(stateVec(:,1));
% cannot find which paper this is from
% if f <= 0.97
%     g = @(f)1.058;
% else
%     g = @(f)1.058 - (0.0155*(f-0.97))/(1.02-f^(1.4));
% end
% alpha = 1.62*exp(0.066*f/(g(f)-f));
% if f <=0.80 || f>=0.995
%     warning('f might be out of range for r_eq computation\n')
% end
% beta_f = exp(0.00077*f/(1.009-f));
% rd = (3*m_s/(4*pi*2165))^(1/3);% Andreas 2005 eq 3.5
% r_eq_fmm = alpha*rd^beta_f;

r_eq_fmm = min_r;

%% tau_r
r_at_tau_r = r_eq+(r_0-r_eq)*exp(-1);
i_star = sum(stateVec(:,1)>r_at_tau_r);

if i_star==length(stateVec(:,1))
    tau_r_fmm = NaN;
else
    % fmm => full microphysical model
    tau_r_fmm = (r_at_tau_r - stateVec(i_star,1))/((stateVec(i_star+1,1) - stateVec(i_star,1))/(time_vec(i_star+1,1)-time_vec(i_star,1)))+time_vec(i_star,1);
end

%% What is used in plotAndreas05_Fig11
T_eq_C = Kelvin2Celsius(T_eq);
T_at_tau_T = T_eq_C+(T_s_0-T_eq_C)*exp(-1);
%% Another possible claculation that will change tau_T_fmm
%     T_at_tau_T = T_eq_fmm+(T_s_0-T_eq_fmm)*exp(-1);
i_star_T = sum(stateVec(:,2)>T_at_tau_T);
% fmm => full microphysical model
tau_T_fmm = (T_at_tau_T - stateVec(i_star_T,2))/((stateVec(i_star_T+1,2) - stateVec(i_star_T,2))/(time_vec(i_star_T+1)-time_vec(i_star_T)))+time_vec(i_star_T);


