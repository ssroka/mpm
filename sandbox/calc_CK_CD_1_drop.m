function [CK, CD] =calc_CK_CD_1_drop(resultsFiles)

%% ------------------------------------------------------
%            Load in the radius and temperature time histories

load(resultsFiles)

%% ----------------------------- setup paths --------------------------
[pathstr,name,ext] = fileparts(mfilename('fullpath')) ;
added_path_flag = 0;
path_to_fxns = which('rho_s.m');
if isempty(path_to_fxns)
fprintf('Using files in %s\n',path_to_fxns(1:end-7))
d = dir(pathstr);
for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        addpath(sprintf('%s/%s',pathstr,d(ii).name))
    end
end
addpath(pathstr)
added_path_flag = 1;
end


helpingAnonFxns;

global Nayar_flag
Nayar_flag = ic.Nayar_flag;

%% ------------------------------------------------------
%           Calculate the sensible and latent heat
%           transferred from the drop

%------------------------------------------------------
% Calculate Volume Flux
%------------------------------------------------------

% if strcmp(ic.inputFile,'Andreas05_Fig11_in')
%     dFdr0 = 1/r_0*r_0; % 1/ (m^2 s^1 ) units
%     Vol_flux = (4.*pi.*ic.r_0.^3/3)*dFdr0;
% elseif strcmp(ic.inputFile,'Andreas08_Fig2_in')
%     % dFdr0 = 3.8*1e-6*U10^(3.4)*5.0*1e-6% Andreas 08
%     dFdr0  = SSGF_Fairall94(ic.r_0*1e6);% argument comes in as micrometers
%     W_u = 3.8*1e-6*ic.U10^(3.4); % whitecap fraction
%     Vol_flux = W_u*dFdr0;
% elseif strcmp(ic.inputfile,'OrtizSuslow')
%     dFdr0  = SSGF_OrtizSuslow(ic.r_0*1e6);% argument comes in as micrometers
%     Vol_flux = 4/3*pi*ic.r_0^3*dFdr0; % units are m^2/(s * micrometers)
    
    Vol_flux =4/3*pi*ic.r_0^3;

% end

%------------------------------------------------------
% reference Andreas 2008 equation 2.4
%------------------------------------------------------
i_star = sum(time_vec<tau_f);
if i_star==length(time_vec)
    T_at_tau_f = T_eq_fmm;
    r_at_tau_f = r_eq_fmm;
else
    T_at_tau_f = (T_s_t(i_star+1)-T_s_t(i_star))/(time_vec(i_star+1)-time_vec(i_star))*(tau_f - time_vec(i_star))+T_s_t(i_star);
    r_at_tau_f = (r_t(i_star+1)-r_t(i_star))/(time_vec(i_star+1)-time_vec(i_star))*(tau_f - time_vec(i_star))+r_t(i_star);
end

Q_L = rho_s(T_s_t(1),r_t(1),ic.m_s,ic.S0/1000,ic.p0).*L_v(T_s_t(1),ic.S0).*(1-(r_at_tau_f/r_t(1))^3).*Vol_flux;
Q_s = rho_s(T_s_t(1),r_t(1),ic.m_s,ic.S0/1000,ic.p0).*c_ps(T_s_t(1),ic.S0,ic.p0).*(T_s_t(1)-T_eq_fmm).*(1-exp(-tau_f/tau_T_fmm)).*Vol_flux;
%%



load('microphysicalConstants.mat','R')

%% ------------------------------------------------------
% 	     Calculate CK
%

% Zweers 2015 equation 7
% H_k = rho_a CK U_L (k_0 - k_L)
%
% H_k    = enthalpy flux
% rho_a  = air density
% CK     = enthalpy flux coefficient
% U_L    = wind speed at height L above the quiescent surface
%        we will generally use 10m
% k      = specific enthalpy (k = k(z))
%        = c_p theta + L_v q
% c_pm   = the specific heat of moist air at constant pressure
% theta  = potential temperature of air
% L_v    = the latent heat of vaporization
% q      = the specific humidity

Hk = Q_s - Q_L;
P0 = 100000; % Pa

T_s_Kelvin = Celsius2Kelvin(ic.T_s);
T_a_Kelvin = Celsius2Kelvin(ic.T_a);

[~, s] = compute_m_w(T_s_t(length(T_s_t)),r_t(length(r_t)),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);

%compute potential temperature
theta_0 = T_s_Kelvin*(P0/ic.p0)^(R/c_pm(ic.T_s,ic.p0,ic.RH));
theta_10 = T_a_Kelvin*(P0/ic.p0)^(R/c_pm(ic.T_a,ic.p0,ic.RH));


% compute specific humidity
% compute mixing ratio
e_0  = ic.RH/100*mb2Pa(e_sat(ic.T_s,Pa2mb(ic.p0)));% vapor pressure
r_v_0 = 0.622*e_0/(ic.p0-e_0); % AMS Glossery website

e_10  = ic.RH/100*mb2Pa(e_sat(ic.T_a,Pa2mb(ic.p0)));% vapor pressure
r_v_10 = 0.622*e_10/(ic.p0-e_10); % AMS Glossery website

q_0 = r_v_0/(1+r_v_0); % AMS Glossery website
q_10 = r_v_10/(1+r_v_10); % AMS Glossery website

k_0 = c_pm(ic.T_s,ic.p0,ic.RH)* theta_0 + L_v(ic.T_s,s*1000)*q_0;
k_10 = c_pm(ic.T_a,ic.p0,ic.RH)* theta_10 + L_v(ic.T_a,s*1000)*q_10;

A_k = 4*pi*ic.r_0^2;
T_k = tau_T_fmm;

CK = Hk/(A_k*T_k*rho_a(Celsius2Kelvin(ic.T_a),ic.p0)*ic.U10*(k_0 - k_10));


%% ------------------------------------------------------
% 	     Calculate CD
%
s_t = zeros(size(T_s_t));
for jj = 1:length(T_s_t)
[~, s_tmp] = compute_m_w(T_s_t(jj),r_t(jj),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);
s_t(jj) = s_tmp;
end
options = odeset('AbsTol',1e-9,'RelTol',1e-9);
rs = rho_s(ic.T_s_0,ic.r_0,ic.m_s,ic.S0/1000,ic.p0);
ra = rho_a(Celsius2Kelvin(ic.T_a),ic.p0);
if ic.r_0 < 100*1e-6
    t_max = 0.5; % seconds
elseif ic.r_0 < 500*1e-6
    t_max = 3; % seconds
elseif ic.r_0 < 1000*1e-6
    t_max = 10; % seconds
else
    t_max = 80; % seconds
end
u0 = 0;

[acel_time_vec, u_vec] = ode45(@(t,u) compute_dudt(t,u,ra,rs,ic.U10,ic.r_0,ic.T_a),[0 t_max],u0,options);
tau_tot = acel_time_vec(max(find(u_vec<ic.U10-ic.U10*0.01)));
tau_f = compute_tauf(ic.U10,ic.T_s_0,ic.r_0,ic.m_s,ic.S0/1000,ic.p0,ic.T_a,ic.maxEr_uf,ic.maxIt);
CD = (rs*4/3*ic.r_0*trapz(acel_time_vec(1:end-1),diff(u_vec)./diff(acel_time_vec)))/(ra*ic.U10^2*min(tau_tot,tau_f));

        
%% ----------------------------- remove user-added paths --------------------------
if added_path_flag
    rm_ASF_paths
end


















