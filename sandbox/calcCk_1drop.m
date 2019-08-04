%% This script will calculate the enthalpy flux coefficient for 1 drop
clear;clc;close all

%% User Input

% define initial microphycial constants and parameters
inputFile = 'Andreas08_Fig2_in'; % Andreas05_Fig11_in Andreas08_Fig2_in Andreas90_Fig09_in CommitteeMeeting_500micormeterdrop

% print the error between the calculated microphysical quantity and the
% quantity stated in the corresponding publication
print_flag = false;

%% ----------------------------- Nayar Functions --------------------------
global Nayar_flag ;
Nayar_flag = false;
if Nayar_flag
	addpath('/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/NayarFxns')
end
%% ----------------------------- BEGIN ------------------------------------

%%  initialize air and sea meterological quantities

eval(inputFile)
fprintf('\nInitial Conditions Set According to: %s \n',inputFile)
fprintf('%s\n\n',(repmat('-',1,38+length(inputFile))))

if exist('microphysicalConstants.mat','file')
    fprintf('\n Microphysical Constants Set\n')
else
    error('\n Please Set Microphysical Constants\n')
end
    fprintf('%s\n\n',(repmat('-',1,38+length(inputFile))))

%% loop through initial radii

r_fail = false(length(r_0_vec),1);
for r_0_ind = 1:length(r_0_vec)
    clearvars -except inputFile print_flag r_fail r_0_ind Nayar_flag
    eval(inputFile)
    
    r_0 = r_0_vec(r_0_ind);
    
    fprintf('\n r_0 = %3.2f %sm %s \n',r_0*1e6,char(181),(repmat('-',1,38+length(inputFile))))
    
    
    
    % -------------------------------------------------------
    %             Sea Spray Generating Function
    % -------------------------------------------------------
    if strcmp(inputFile,'Andreas05_Fig11_in')
        dFdr0 = 1/r_0*r_0; % 1/ (m^2 s^1 ) units
        Vol_flux = (4.*pi.*r_0.^3/3)*dFdr0;
    elseif strcmp(inputFile,'Andreas08_Fig2_in')
        % dFdr0 = 3.8*1e-6*U10^(3.4)*5.0*1e-6% Andreas 08
        dFdr0  = SSGF_Fairall94(r_0*1e6);% argument comes in as micrometers
        dFdr0_vec(r_0_ind)=dFdr0;
        W_u = 3.8*1e-6*U10^(3.4); % whitecap fraction
        Vol_flux = W_u*dFdr0;
    elseif strcmp(inputfile,'Exp1')
        dFdr0  = SSGF_OrtizSuslow(r_0);% argument comes in as micrometers
        dFdr0_vec(r_0_ind)=dFdr0;
        Vol_flux = 4/3*pi*r_0^3*dFdr0; % units are m^2/(s * micrometers) 

    end
    
    % -------------------------------------------------------
    %             Calculate Micro-physical Endpoints
    % -------------------------------------------------------
    helpingAnonFxns;
    s0    = S0/1000; % kg/kg
    T_s   = T_s_0;  % deg Celsius
    [m_s] = compute_initial_drop_conditions(T_s_0,S0,r_0,p0);
    
    storeICs;
    
    compute_Teq;
    
    compute_req;
    
    compute_tauT;
    
    compute_taur;
    
    t_final = 1050;
    
    integrate_r_T;
    
    collect_mm_values;
    
    compute_tauf;
    
    print_comparisons;
    
    saveScript;
end
%{
return
    % --------------------------------------------------------
    %             Evolve r and T
    % --------------------------------------------------------
    timesteps = [0,logspace(-2,1,200)]'; % see Andreas 2005, Fig 1
    
    
    % --------------------------------------------------------
    %                          Qs
    % --------------------------------------------------------
    load('microphysicalConstants.mat','c_ps')
    T_eq_fmm = min(T_s_t);
    T_eq_C = T_eq_fmm;
    
    if ~exist('tau_T_fmm','var')
        T_at_tau_T = T_eq_C+(T_s_0-T_eq_C)*exp(-1);
        i_star_T = sum(stateVec(:,2)>T_at_tau_T);
        % fmm => full microphysical model
        tau_T_fmm = (T_at_tau_T - stateVec(i_star_T,2))/((stateVec(i_star_T+1,2) - stateVec(i_star_T,2))/(time_vec(i_star_T+1)-time_vec(i_star_T)))+time_vec(i_star_T);
    end
    if ~exist('r_tau_f','var')
        i_star_f = sum(time_vec<tau_f);
        % fmm => full microphysical model
        r_tau_f = (r_t(i_star_f+1)-r_t(i_star_f))/(time_vec(i_star_f+1)-time_vec(i_star_f))*(time_vec(i_star_f)-tau_f)+r_t(i_star_f);
    end
    
    
    s_t = m_s./(m_s+compute_m_w(T_s_t,r_t,m_s));
    % TODO: check rho_s, need to actualy interpolate
    Q_s = rho_s(T_s_0,r_0,m_s,s0).*c_ps.*(T_s_0-T_eq_C).*(1-exp(-tau_f/tau_T_fmm)).*Vol_flux;
    Q_s_r_0(r_0_ind) = Q_s(1);
    
    
    % --------------------------------------------------------
    %                          QL
    % --------------------------------------------------------
    % radius at t= tau_f =
    r_f = exp(-tau_f/tau_r)*(r_0-r_eq)+r_eq;
    % TODO: check rho_s, need to actualy interpolate
    
    Q_L = rho_s(T_s_0,r_0,m_s,s0).*L_v(T_s_0).*(1-(r_f/r_0).^3)*Vol_flux;
    
    Q_L_r_0(r_0_ind) = Q_L(1);
    
    % --------------------------------------------------------
    %                      Estimate Ck
    % --------------------------------------------------------
    Hk_t = (Q_s(1)+Q_L(1));
    % Calculate H_k
    % specific humidity
    % q = 0.622*RH*rho_s/(rho_a - rho_s);
    epsilon = 0.622;
    % q = 1000*(epsilon*(e_sat_Pa)/((p0)-((1-epsilon)*(e_sat_Pa))));
    mixing_ratio = epsilon*e_sat_Pa/(p0-e_sat_Pa);
    q = mixing_ratio/(mixing_ratio+1);
    
    % just a linear function that agrees with the aformentioned temperatures assuming the pressure is always the reference pressure
    theta = @(z)(T_a-T_s)/10*z+T_s;
    theta_0_K = Celsius2Kelvin(theta(0));
    theta_10_K = Celsius2Kelvin(theta(10));
    
    k_0_t = c_p*theta_0_K-L_v_t*q;
    k_10_t = c_p*theta_10_K-L_v_t*q;
    
    C_k = Hk_t./(rho_a*U10*(k_0_t-k_10_t));
    % print CK
    % fprintf('max C_K = %g\nmin C_K = %g\n',max(C_k),min(C_k))
    
    % --------------------------------------------------------
    %                      Estimate CD
    % --------------------------------------------------------
    % 5/3 comes from added mass
    % C_D = ((4/3)*rho_s*r_t*U10*(5/3))/((p0+0.5*rho_a*U10^2)*2*rho_a*U10^2)
    
    % C_D = 4/3*pi*r_t.^3.*rho_s_t/(rho_a*U10);
    C_D = (4/3+2/3)*pi*r_t.^3.*rho_s_t/(rho_a*U10);% +2/3 to acount for the added mass
    
    % print CD
    fprintf('max C_D = %g\nmin C_D = %g\n',max(C_D),min(C_D))
    
    
    % scale factor
    scale = max(C_D(2))/max(C_k(2));
    C_k_scale = (scale*Hk_t)./(rho_a*U10*(k_0_t-k_10_t));
    
    % plot ratio Ck to CD
    plot_CK_CD_ratio;
    
    plot_CK_CD_unscaled;
    
    plot_CK_CD_scaled;
    
    
end
%}
%% plotting

% plot_Q_S_Q_L;
% figure(1)
% plot(r_0_vec,y_t,'displayname',num2str(r_0))
% legend('-dynamiclegend')
% xlabel('r_0')
% ylabel('y')
% hold on
% figure(2)
% loglog(r_0_vec,Q_s_r_0,'r','displayname','Q_s')
% hold on
% loglog(r_0_vec,Q_L_r_0,'bs','displayname','Q_L')
% legend('-dynamiclegend')

