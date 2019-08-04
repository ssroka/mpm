function [] =study_indep_vars_real_data(resultsFiles,save_new_data_flag,U10,next_r,DT)

%% ------------------------------------------------------
%            Load in the radius and temperature time histories

load(resultsFiles)

%% ----------------------------- setup paths --------------------------
[pathstr,name,ext] = fileparts(mfilename('fullpath')) ;
added_path_flag = 0;
path_to_fxns = which('rho_s.m');
current_dir = pwd;
if isempty(path_to_fxns)
    loc_of_sandbox = current_dir(1:strfind(current_dir,'sandbox')+6);
    %     cd(loc_of_sandbox)
    d = dir(loc_of_sandbox);
    dirs_to_add = {d.name};
    for ii = 1:length(dirs_to_add)
        % do not include hidden directories
        if ~strcmp(dirs_to_add{ii}(1),'.') && isdir(sprintf('%s/%s',loc_of_sandbox,dirs_to_add{ii}))
            addpath(sprintf('%s/%s',loc_of_sandbox,dirs_to_add{ii}))
        end
    end
    addpath(loc_of_sandbox)
    added_path_flag = 1;
else
    fprintf('Using files in %s\n',path_to_fxns(1:end-7))
end


helpingAnonFxns;
load('microphysicalConstants.mat','R')

global Nayar_flag
Nayar_flag = ic.Nayar_flag;

%% ------------------------------------------------------
%           Calculate the sensible and latent heat
%           transferred from the drop

[~,n,~]=fileparts(resultsFiles);
save_dat_file_test = strcat('./',n,'_energyX_',num2str(U10),'_',num2str(DT),'.mat');
%% ----------------------------- Begin --------------------------
if save_new_data_flag
    if exist(save_dat_file_test, 'file') == 2
        save_dat_file = sprintf('%s_%s',n,datestr(now,'yyyymmdd_HHMMSS'));
        clc;
        fprintf('data file %s already exists,\nchanging name to %s \n',save_dat_file_test,save_dat_file)
    else
        save_dat_file=save_dat_file_test;
    end
    
    Q_s = zeros(length(time_vec)-1,1);
    Q_L = zeros(length(time_vec)-1,1);
    C_D = zeros(length(time_vec)-1,1);
    C_K = zeros(length(time_vec)-1,1);
    k_0 = zeros(length(time_vec)-1,1);
    k_10 = zeros(length(time_vec)-1,1);
    C_K_int_denom = zeros(length(time_vec)-1,1);
    m_w = zeros(2,1);
    s   = zeros(2,1);
    L_v = zeros(2,1);
    Volume = zeros(2,1);
    
    g = 9.81; % m s^-2
    T_a_Kelvin = Celsius2Kelvin(ic.T_a);
    
    %compute potential temperature
    theta_10 = Celsius2Kelvin(ic.T_a)*(ic.p0/ic.p0)^(R/c_pm(ic.T_a,ic.p0,ic.RH));
    e_10  = ic.RH/100*mb2Pa(e_sat(ic.T_a,Pa2mb(ic.p0)));% vapor pressure
    r_v_10 = 0.622*e_10/(ic.p0-e_10); % AMS Glossery website
    q_10 = r_v_10/(1+r_v_10); % AMS Glossery website
    
    fprintf('starting : %s\n',n)
    toast_points = round(length(time_vec)*(1:3)/4);
    tic
    
    % compute speed of drop
    t_max = max(time_vec); % seconds
    u0    = 0; % m/s starting from rest
    ra = rho_a(Celsius2Kelvin(ic.T_a),ic.p0);
    rs = rho_s(ic.T_s_0,ic.r_0,ic.m_s,ic.S0/1000,ic.p0);
    options = odeset('AbsTol',1e-9,'RelTol',1e-9);
    [acel_time_vec, u_vec_tmp] = ode45(@(t,u) compute_dudt(t,u,ra,rs,U10,ic.r_0,ic.T_a),[0 t_max],u0,options);
    u_vec = interp1(acel_time_vec,u_vec_tmp,time_vec);
    
    
    
    for t = 1:length(time_vec)-1
        if ismember(t,toast_points)
            fprintf('%d%% done (%f seconds) \n',round(t/length(time_vec)*100),toc)
            tic
        end

        if next_r
            % Step 1: compute average saltiness [kg/kg]
            m_w(1) = compute_m_w(T_s_t(t),r_t(t),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);
            s(1) = ic.m_s/(m_w(1)+ic.m_s);
            m_w(2) = compute_m_w(T_s_t(t+1),r_t(t+1),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);
            s(2) = ic.m_s/(m_w(2)+ic.m_s);
            
            % Step 2: compute average volume [m^3]
            Volume(1) = 4/3*pi*r_t(t).^3;
            Volume(2) = 4/3*pi*r_t(t+1).^3;
            Volume_avg = mean(Volume);
            
            % Step 3: compute specific heat [J kg^-1 K^-1]
            cps_avg = 0.5*(SW_SpcHeat(T_s_t(t+1),'C',s(2),'w',ic.p0,'Pa')+SW_SpcHeat(T_s_t(t),'C',s(1),'w',ic.p0,'Pa'));
            cp_avg = 0.5*(c_pm(T_s_t(t+1),ic.p0, ic.RH)+c_pm(T_s_t(t),ic.p0, ic.RH));
            
            % Step 4: compute average density over this time [kg m^-3]
            rho_s_avg = 0.5*(rho_s(T_s_t(t+1),r_t(t+1),ic.m_s,s(2),ic.p0) + rho_s(T_s_t(t),r_t(t),ic.m_s,s(1),ic.p0));
            
            % Step 5: compute Q_s
            Q_s(t) = rho_s_avg*Volume_avg*cps_avg*(Celsius2Kelvin(T_s_t(1))-Celsius2Kelvin(T_s_t(t)));
            
            % Step 6: compute average L_v [J/kg]
            L_v(1) = SW_LatentHeat(T_s_t(t),'C',s(1),'w');
            L_v(2) = SW_LatentHeat(T_s_t(t+1),'C',s(2),'w');
            L_v_avg = mean(L_v);
            
            % Step 7: compute Q_L
            Q_L(t) = rho_s_avg*4/3*pi*r_t(1)^3*(1-(r_t(t+1)/r_t(1))^3)*L_v_avg;
            
            
            T_s_Kelvin = Celsius2Kelvin(T_s_t(t));
            theta_0 = T_s_Kelvin*(ic.p0/ic.p0)^(R/c_pm(T_s_t(t),ic.p0,100));
            e_0  = 100/100*mb2Pa(e_sat(T_s_t(t),Pa2mb(ic.p0)));% vapor pressure
            r_v_0 = 0.622*e_0/(ic.p0-e_0); % AMS Glossery website
            q_0 = r_v_0/(1+r_v_0); % AMS Glossery website
            
            %         T_s_Kelvin_init = Celsius2Kelvin(ic.T_s);
            %         theta_0 = T_s_Kelvin_init*(ic.p0/ic.p0)^(R/c_pm(ic.T_s,ic.p0,100));
            %         e_0  = 100/100*mb2Pa(e_sat(ic.T_s,Pa2mb(ic.p0)));% vapor pressure
            %         r_v_0 = 0.622*e_0/(ic.p0-e_0); % AMS Glossery website
            %         q_0 = r_v_0/(1+r_v_0); % AMS Glossery website
            
            
            k_0(t)  = cp_avg* theta_0  + L_v_avg*q_0;
            k_10(t) = cp_avg* theta_10 + L_v_avg*q_10;
            
            C_K(t) = (Q_s(t)-Q_L(t));
            C_K_int_denom(t) = (4*pi*r_t(t)^2*rho_a(T_a_Kelvin,ic.p0)*(k_10(t)-k_0(t)));
        else
            load(save_dat_file_test,'C_K')
        
        end
        
        h = u_vec(t).^2/g/2;
        C_D(t) = rho_s_avg*4/3*pi*r_t(t).^3*(g*h+ 0.5*u_vec(t).^2);
        
    end
    C_K_int_num = (cumsum(Q_s.*diff(time_vec)./time_vec(2:end))+cumsum(Q_L.*diff(time_vec)./time_vec(2:end)));
    
    
    %% save data file
    
    save(save_dat_file)
    
else
    save_dat_file = save_dat_file_test;
    load(save_dat_file)
    
end



%%


t_lims = find(time_vec<20);
t_lims_p1 = [t_lims;t_lims(end)+1];
t_lims_m1 = [t_lims_p1(2:end)];


not_yet_evap = Q_s-Q_L>0;

% figure(2)
% subplot(3,1,1)
% semilogx(time_vec(not_yet_evap), C_K(not_yet_evap),'linewidth',3)
% hold on
% title('Net Thermal Energy Given to Atmosphere [J]')
% subplot(3,1,2)
% semilogx(time_vec(not_yet_evap), C_D(not_yet_evap),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m, U = %2d m/s',round(ic.r_0*1e6),char(956),U10))
% hold on
% lh = legend('-dynamiclegend');
% set(lh,'location','west');
% title('Net Mechanical Energy Taken from Atmosphere [J]')
% subplot(3,1,3)
% semilogx(time_vec(not_yet_evap), C_K(not_yet_evap)./C_D(not_yet_evap),'linewidth',3)
% title('Net Thermal Energy Given / Net Mechanical Energy Taken')
% xlabel('t [s]')
% hold on
% set(gcf,'color','w','position',[440   253   497   545])
% text(2e-3,3,'T_a = 18 C')
% text(2e-3,2.5,'T_s = 20 C')
% text(2e-3,2,'U_{10} = 50 m/s')
% text(2e-3,1.5,'RH = 90 %')

figure(2)
subplot(3,1,1)
plot3(ic.r_0*ones(sum(not_yet_evap),1),time_vec(not_yet_evap), C_K(not_yet_evap),'linewidth',3)
hold on
title('Net Thermal Energy Given to Atmosphere [J]')
subplot(3,1,2)
plot3(ic.r_0*ones(sum(not_yet_evap),1),time_vec(not_yet_evap), C_D(not_yet_evap),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m, U = %2d m/s, \Delta T = %1.1f',round(ic.r_0*1e6),char(956),U10,DT))
hold on
lh = legend('-dynamiclegend');
set(lh,'location','west');
title('Net Mechanical Energy Taken from Atmosphere [J]')
subplot(3,1,3)
plot3(ic.r_0*ones(sum(not_yet_evap),1),time_vec(not_yet_evap), C_K(not_yet_evap)./C_D(not_yet_evap),'linewidth',3)
title('Net Thermal Energy Given / Net Mechanical Energy Taken')
ylabel('t [s]')
xlabel('r_0 [\mu m]')
hold on
set(gcf,'color','w','position',[440   253   497   545])
text(2e-3,3,'T_a = 18 C')
text(2e-3,2.5,'T_s = 20 C')
text(2e-3,2,'U_{10} = 50 m/s')
text(2e-3,1.5,'RH = 90 %')


if added_path_flag
    rm_ASF_paths
end
%{

addpath ~/Documents/MATLAB/util/export_fig/
export_fig(gcf,'QsQL.pdf','-m2')
export_fig(gcf,'QsQL.fig')
copyfile('QsQL.pdf','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')
copyfile('QsQL.fig','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')


addpath ~/Documents/MATLAB/util/export_fig/
export_fig(gcf,'Qs_minus_QL.pdf','-m2')
export_fig(gcf,'Qs_minus_QL.fig')
copyfile('Qs_minus_QL.pdf','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')
copyfile('Qs_minus_QL.fig','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')


addpath ~/Documents/MATLAB/util/export_fig/
export_fig(gcf,'CKCD_new.pdf','-m2')
export_fig(gcf,'CKCD_new.fig','-m2')
copyfile('CKCD_new.pdf','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')
copyfile('CKCD_new.fig','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')
copyfile('CKCD_new.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')
copyfile('CKCD_new.fig','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')
    
addpath ~/Documents/MATLAB/util/export_fig/
export_fig(gcf,'max_CKCD_new.pdf','-m2')
export_fig(gcf,'max_CKCD_new.fig','-m2')
copyfile('max_CKCD_new.pdf','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')
copyfile('max_CKCD_new.fig','~/MIT/Research/EmanuelGroup/CommitteeMeeting2/')
copyfile('max_CKCD_new.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')
copyfile('max_CKCD_new.fig','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')
    
    




%}















