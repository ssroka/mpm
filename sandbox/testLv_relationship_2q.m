

for RH = 1:11
T_s_Kelvin = Celsius2Kelvin(T_s_t(t));
theta_10 = T_s_Kelvin*(ic.p0/ic.p0)^(R/c_pm(T_s_t(t),ic.p0,100));
e_10  = (89+RH)/100*mb2Pa(e_sat(T_s_t(t),Pa2mb(ic.p0)));% vapor pressure
r_v_10 = 0.622*e_10/(ic.p0-e_10); % AMS Glossery website
q_10(RH) = r_v_10/(1+r_v_10); % AMS Glossery website
end
plot(90:100,q_10)









        