function [dsdt] = compute_dsdt(t,stateVec,RH,m_s,T_a,P_mb,t_final,R,M_w,maxEr_s,maxIt)
helpingAnonFxns;
P_Pa = mb2Pa(P_mb);

progress_fxn(t/t_final);

% full microphysical model evolving both r and T together

% T_a, the an=mbient atmospheric temperature in Celsius, is assumed to not
% change
% state vector stateVec = [r ; T]'

r   = stateVec(1);
T_s = stateVec(2);

T_a_Kelvin = Celsius2Kelvin(T_a);
T_s_Kelvin = Celsius2Kelvin(T_s);

% drdt = equation 2.1 in Andreas 2005

% -------------------------------------------------------------------
%    		  		f
% -------------------------------------------------------------------
f = RH/100;

% -------------------------------------------------------------------
%    		  		m_w
% -------------------------------------------------------------------
m_w = compute_m_w(T_s,r,m_s,P_Pa,maxEr_s,maxIt);

% -------------------------------------------------------------------
%    		  		s
% -------------------------------------------------------------------
s = m_s/(m_w+m_s);

% -------------------------------------------------------------------
%    		  		e_sat
% -------------------------------------------------------------------
e_sat_Pa =  mb2Pa(e_sat(T_a,P_mb));

% -------------------------------------------------------------------
%    		  		drdt
% -------------------------------------------------------------------

% constants: 		f, R, M_w, T_a_Kelvin, e_sat_Pa
% vary in time: 	y,r,rho_s,Phi_s, r, rho_s, D_w_prime, k_a_prime(r(t)), L_v(T_s(t))

drdt   =(((f-1)-y(T_s,T_a,r,m_s,s,P_Pa))./r)/(rho_s(T_s,r,m_s,s,P_Pa).*R.*T_a_Kelvin./(D_w_prime(T_s,P_mb,r).*M_w.*e_sat_Pa)+rho_s(T_s,r,m_s,s,P_Pa).*L_v(T_s,s*1000)./(k_a_prime(T_s,P_mb,r).*T_a_Kelvin).*(L_v(T_s,s*1000).*M_w./(R.*T_a_Kelvin)-1));


e_sat_T_a_Pa= mb2Pa(e_sat(T_a,P_mb));

e_sat_T_s_Pa= mb2Pa(e_sat(T_s,P_mb));

rho_v = f*M_w*e_sat_T_a_Pa/(R*T_a_Kelvin);
rho_vr = M_w*e_sat_T_s_Pa*exp(y(T_s,T_a,r,m_s,s,P_Pa))/(R*T_s_Kelvin);

% -------------------------------------------------------------------
%    		  		dTdt
% -------------------------------------------------------------------

dTdt = 3./(rho_s(T_s,r,m_s,s,P_Pa)*c_ps(T_s,s*1000,P_Pa)*r^2)*(k_a_prime(T_s,P_mb,r)*(T_a_Kelvin - T_s_Kelvin)+L_v(T_s,s*1000)*D_w_prime(T_s,P_mb,r)*(rho_v - rho_vr));


% -------------------------------------------------------------------
%    		  		dsdt
% -------------------------------------------------------------------

dsdt = [drdt;dTdt];
