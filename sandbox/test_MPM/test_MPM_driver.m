% this is a driver file that will test the MPM code against the example in
% Andreas 2005 Figure 11 to confirm that the integration yeilds the same
% result

% if the way the driver calls the solver script dropletEvolution, or the
% obligations of the driver to define variables/add paths/ etc. change then
% those updates must be made manually to this file

% The directory containing this driver should be at the same level 
% as inputFiles

clear ;clc;close all

if exist('./results/Test_MPM','dir')
    error('Please remove the current test directory \n %s/results/Test_MPM',pwd)
end

%% User Input

% define initial microphycial constants and parameters
inputFile = 'test_MPM_in'; 
%% ----------------------------- Nayar Functions --------------------------
global Nayar_flag ;

%% ----------------------------- setup paths --------------------------
pwd_tmp = pwd;
dirstr = pwd_tmp(1:strfind(pwd,'sandbox')+length('sandbox')-1);
d = dir(dirstr);
addpath(dirstr)
for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        addpath(sprintf('%s/%s',dirstr,d(ii).name))
    end
end

%% ----------------------------- Run Microphysical Model --------------------------
dropletEvolution;

%% ----------------------------- remove user-added paths --------------------------

% rm_ASF_paths;

%% ----------------------------- load solutions --------------------------
% clear;clc
vars_2_compare = {'T_eq_fmm','r_eq_fmm','tau_T_fmm','tau_r_fmm','T_eq',...
                  'r_eq','tau_T','tau_r'};
              
load('./results/Test_MPM/r_0_100_Test_MPM.mat',vars_2_compare{:})

for jj = 1:length(vars_2_compare)
    eval(sprintf('%s_new = %s;',vars_2_compare{jj},vars_2_compare{jj}));
    eval(sprintf('clear %s;',vars_2_compare{jj}));
end

load('Andreas05_test_MPM.mat')
T_s_t = stateVec(:,2);
T_a = Kelvin2Celsius(T_a_Kelvin);
T_eq_vec = zeros(size(T_s_t));
rho_v  = (M_w*e_sat_Pa/(R*T_a_Kelvin))*f;
rho_vr = (M_w*mb2Pa(e_sat(T_s,P_mb))/(R*Celsius2Kelvin(T_s)))*exp(y(T_s,T_a,r_0,m_s,s(r_0),ic.p0));

T_eq_vec(1) = T_a+L_v(T_s)*D_w_prime(T_s,P_mb,r_0)*(rho_v-rho_vr)/k_a_prime(T_s,P_mb,r_0);
e_sat_Pa = mb2Pa(e_sat(T_a,P_mb));
rho_v  = (M_w*e_sat_Pa/(R*T_a_Kelvin))*f;
for i = 1:length(T_s_t)-1
    rho_vr = (M_w*mb2Pa(e_sat(T_eq_vec(i),P_mb))/(R*Celsius2Kelvin(T_eq_vec(i))))*exp(y(T_eq_vec(i),T_a,r_t(i+1),m_s,s(r_t(i+1)),ic.p0));
    T_eq_vec(i+1) = 0.5*(T_eq_vec(i) + T_a + L_v(T_s_t(i))*D_w_prime(T_s_t(i),P_mb,r_t(i))*(rho_v-rho_vr)/k_a_prime(T_s_t(i),P_mb,r_t(i)));
end
kk = find(time_vec<tau_T);
T_eq_fmm = T_eq_vec(max(kk)+1);

r_eq_fmm = min(stateVec(:,1));


% test microphysical model quantities
fmm_vars_match = abs(abs([T_eq_fmm_new;r_eq_fmm_new;tau_T_fmm_new;tau_r_fmm_new]-[T_eq_fmm;r_eq_fmm;tau_T_fmm;tau_r_fmm])./[T_eq_fmm;r_eq_fmm;tau_T_fmm;tau_r_fmm])>1e-6;
if any(fmm_vars_match)
    fprintf('Failed: Four endpoints from full model\n')
    fprintf('Failed in variable %s\n',vars_2_compare{[fmm_vars_match;false(4,1)]})
else
    fprintf('Passed: Four endpoints from full model\n')
end

% test microphysical model quantities
af_vars_match = abs(abs([T_eq_new;r_eq_new;tau_T_new;tau_r_new]-[T_eq;r_eq;tau_T;tau_r])./[T_eq;r_eq;tau_T;tau_r])>1e-6;
if any(af_vars_match)
    fprintf('Failed: Four endpoints from approximations\n')
    fprintf('Failed in variable %s\n',vars_2_compare{[false(4,1) af_vars_match]})
else
    fprintf('Passed: Four endpoints from approximations\n')
end


