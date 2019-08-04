% use ode45 to integrate drdt to arrive at r
% then fit a high order polynomial to this r and
% and take two derivatives in time to solve for d2rdt2 to use in solving
% for tau_r
options = odeset('AbsTol',1e-15,'RelTol',2.22045e-14);
T_eq_C = Kelvin2Celsius(T_eq);
P_mb = Pa2mb(p0);
t_final_r = t_final;
while ~exist('r_vec','var')
    try
        [time_vec_r, r_vec] = ode45(@(t,r) compute_drdt(t,r,tau_T,T_eq_C,T_s_0,T_a,RH,P_mb,m_s,maxEr_s,maxIt),[0 t_final_r],r_0,options);
    catch
        t_final_r = 0.5*t_final_r;
    end
end






































