function [CK, CD] =new_calc_CK_CD(resultsFiles,save_new_data_flag)

%% ------------------------------------------------------
%            Load in the radius and temperature time histories

load(resultsFiles)

%% ----------------------------- setup paths --------------------------
[pathstr,name,ext] = fileparts(mfilename('fullpath')) ;
added_path_flag = 0;
path_to_fxns = which('rho_s.m');
if isempty(path_to_fxns)
    d = dir(pathstr);
    for ii = 1:length(d)
        % do not include hidden directories
        if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
            addpath(sprintf('%s/%s',pathstr,d(ii).name))
        end
    end
    addpath(pathstr)
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
save_dat_file = strcat('./',n,'_CKCDQsQL_Ts_t_dat.mat');
%% ----------------------------- Begin --------------------------
if save_new_data_flag
    if exist(save_dat_file, 'file') == 2
        save_dat_file = sprintf('%s_CKCDQsQL_dat_%s',n,datestr(now,'yyyymmdd_HHMMSS'));
        clc;
        fprintf('data file %s already exists,\nchanging name to %s \n',strcat(n,'_CKCDQsQL_dat'),save_dat_file)
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
    
    T_a_Kelvin = Celsius2Kelvin(ic.T_a);
    
           %compute potential temperature
        theta_10 = Celsius2Kelvin(ic.T_a)*(ic.p0/ic.p0)^(R/c_pm(ic.T_a,ic.p0,ic.RH));
        e_10  = ic.RH/100*mb2Pa(e_sat(ic.T_a,Pa2mb(ic.p0)));% vapor pressure
        r_v_10 = 0.622*e_10/(ic.p0-e_10); % AMS Glossery website
        q_10 = r_v_10/(1+r_v_10); % AMS Glossery website
    
    fprintf('starting : %s\n',n)
    toast_points = round(length(time_vec)*(1:3)/4);
    tic
    for t = 1:length(time_vec)-1
        if ismember(t,toast_points)
        fprintf('%d%% done (%f seconds) \n',round(t/length(time_vec)*100),toc)
        tic
        end
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
        
        g = 9.81;U10=50;h = U10.^2/g/2;
        C_D(t) = rho_s_avg*4/3*pi*r_t(t).^3*(g*h+ 0.5*U10.^2);
        
    end
    C_K_int_num = (cumsum(Q_s.*diff(time_vec)./time_vec(2:end))+cumsum(Q_L.*diff(time_vec)./time_vec(2:end)));
    
    
    %% calculate CD
    
    
    % for t = 1:length(time_vec)-1
    %     % Step 1: compute average saltiness [kg/kg]
    %     m_w(1) = compute_m_w(T_s_t(t),r_t(t),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);
    %     s(1) = ic.m_s/(m_w(1)+ic.m_s);
    %     m_w(2) = compute_m_w(T_s_t(t+1),r_t(t+1),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);
    %     s(2) = ic.m_s/(m_w(2)+ic.m_s);
    %
    %     % Step 2: compute average volume [m^3]
    %     Volume(1) = 4/3*pi*r_t(t)^3;
    %     Volume(2) = 4/3*pi*r_t(t+1)^3;
    %     Volume_avg = mean(Volume);
    %
    %     % Step 4: compute average density over this time [kg m^-3]
    %     rho_s_avg = 0.5*(rho_s(T_s_t(t+1),r_t(t+1),ic.m_s,s(2),ic.p0) + rho_s(T_s_t(t),r_t(t),ic.m_s,s(1),ic.p0));
    %
    % end
    
    %% save data file
    
    save(save_dat_file)
    
else
    
    load(save_dat_file)
    
end



%%


t_lims = find(time_vec<20);
t_lims_p1 = [t_lims;t_lims(end)+1];
t_lims_m1 = [t_lims_p1(2:end)];


figure(1)
yyaxis left
semilogx(time_vec(t_lims),cumsum(Q_s(t_lims).*diff(time_vec(t_lims_p1))./time_vec(t_lims_m1)),'linewidth',3,'displayname',sprintf('Q_s r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
ylabel('Integrated Q_s [J]')
% set(gca,'ylim',[0 1.2]*1e-4,'xlim',[1e-3 1e2])
yyaxis right
semilogx(time_vec(t_lims),abs(cumsum(Q_L(t_lims).*diff(time_vec(t_lims_p1))./time_vec(t_lims_m1))),'linewidth',3,'displayname',sprintf('Q_L r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
ylabel('Integrated Q_L [J]')
set(gca,'ylim',[0 1.2]*1e-4,'xlim',[1e-3 1e2])
xlabel('t [s]')
set(gcf,'color','w','position',[440   300   737   498])
lh = legend('-dynamiclegend');
set(gca,'fontsize',20)
set(lh,'location','westoutside','fontsize',12)


figure(11)
subplot(121)
yyaxis left
semilogx(time_vec(t_lims),Q_s(t_lims),'linewidth',3,'displayname',sprintf('Q_s r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
ylabel(' Q_s [J]')
set(gca,'ylim',[0 1.2]*1e-4,'xlim',[1e-3 1e2])
yyaxis right
semilogx(time_vec(t_lims),Q_L(t_lims),'linewidth',3,'displayname',sprintf('Q_L r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
ylabel(' Q_L [J]')
set(gca,'ylim',[0 1.2]*1e-4,'xlim',[1e-3 1e2])
xlabel('t [s]')
subplot(122)
semilogx(time_vec(t_lims),(Q_s(t_lims)-Q_L(t_lims)),'linewidth',3,'displayname',sprintf('Q_L r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
set(gcf,'color','w','position',[440   300   737   498])
lh = legend('-dynamiclegend');
set(gca,'fontsize',20)
set(lh,'location','westoutside','fontsize',12)

figure(12)
semilogx(time_vec(t_lims),cumsum(Q_s(t_lims).*diff(time_vec(t_lims_p1))./time_vec(t_lims_m1))+cumsum(Q_L(t_lims).*diff(time_vec(t_lims_p1))./time_vec(t_lims_m1)),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
xlabel('t [s]')
title('Net Heating [J]')
% set(gca,'ylim',[-4 4]*1e-5,'xlim',[1e-3 1e2])
set(gcf,'color','w','position',[ 440   384   616   414])
lh = legend('-dynamiclegend');
set(gca,'fontsize',20)
set(lh,'location','westoutside','fontsize',12)
%

not_yet_evap = Q_s-Q_L>0;

figure(2)
subplot(3,1,1)
semilogx(time_vec(not_yet_evap), C_K(not_yet_evap),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)))
hold on
lh = legend('-dynamiclegend');
set(lh,'location','west');
title('Net Thermal Energy Given to Atmosphere [J]')
subplot(3,1,2)
semilogx(time_vec(not_yet_evap), C_D(not_yet_evap),'linewidth',3)
hold on
title('Net Mechanical Energy Taken from Atmosphere [J]')
subplot(3,1,3)
semilogx(time_vec(not_yet_evap), C_K(not_yet_evap)./C_D(not_yet_evap),'linewidth',3)
title('Net Thermal Energy Given / Net Mechanical Energy Taken')
xlabel('t [s]')
hold on
set(gcf,'color','w','position',[440   253   497   545])
text(2e-3,3,'T_a = 18 C')
text(2e-3,2.5,'T_s = 20 C')
text(2e-3,2,'U_{10} = 50 m/s')
text(2e-3,1.5,'RH = 90 %')


C_K_3 = C_K;
C_D_3 = C_D;
C_K_3(~not_yet_evap) = 0;
LT_300 = time_vec<300;
figure(3)
subplot(3,1,1)
semilogx(time_vec(LT_300), C_K_3(LT_300),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)))
hold on
lh = legend('-dynamiclegend');
set(lh,'location','west');
title('Q_s - Q_L [J]')
subplot(3,1,2)
semilogx(time_vec(LT_300), C_D_3(LT_300),'linewidth',3)
hold on
title('1/2 m U^2 + m g h [J]')
subplot(3,1,3)
semilogx(time_vec(LT_300), C_K_3(LT_300)-C_D_3(LT_300),'linewidth',3)
title('Q_s - Q_L - 1/2 m U^2 - m g h [J]')
xlabel('t [s]')
hold on
set(gcf,'color','w','position',[440   253   497   545])
text(2e-3,3,'T_a = 18 C')
text(2e-3,2.5,'T_s = 20 C')
text(2e-3,2,'U_{10} = 50 m/s')
text(2e-3,1.5,'RH = 90 %')


%  time of max heat transfer
[val,ind] = max(C_K(not_yet_evap)./C_D(not_yet_evap));


fig21_dat = matfile('fig21.mat','Writable',true);
try
    fig21_dat.r = [fig21_dat.r;ic.r_0];
    fig21_dat.CK2CD = [fig21_dat.CK2CD;val];
    fig21_dat.CK2CD_t = [fig21_dat.CK2CD_t;time_vec(ind)];
catch
    fig21_dat.r = ic.r_0;
    fig21_dat.CK2CD = val;
    fig21_dat.CK2CD_t = time_vec(ind);
end
% if length(fig21_dat.r)<1
%     fig21_dat.r = ic.r_0;
%     fig21_dat.CK2CD = max(C_K(not_yet_evap)./C_D(not_yet_evap));
% else
%     fig21_dat.r(end+1) = ic.r_0;
%     fig21_dat.CK2CD(end+1) = max(C_K(not_yet_evap)./C_D(not_yet_evap));
% end
figure(21)
clf
load('fig21.mat')
yyaxis left
plot(fig21_dat.r*1e6,fig21_dat.CK2CD,'*-','linewidth',3)
ylabel(sprintf('max(Thermal Energy Given /Mechanical Energy Taken) [J/J]'))
yyaxis right
plot(fig21_dat.r*1e6,fig21_dat.CK2CD_t,'s-','linewidth',3)
ylabel(sprintf('time at which \nmax(Thermal Energy Given / Mechanical Energy Taken)\n occurs [s]'))
xlabel(sprintf('r_0 %s m',char(956)))
set(gcf,'color','w','position',[ 440   384   616   414])
text(55,3.2,'T_a = 18 C')
text(55,3,'T_s = 20 C')
text(55,2.8,'U_{10} = 50 m/s')
text(55,2.6,'RH = 90 %')


%
% CK_gt_0 = C_K_int_num>0;
%
% figure(3)
% subplot(3,1,1)
% semilogx(time_vec(CK_gt_0), C_K(CK_gt_0),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_K')
% subplot(3,1,2)
% semilogx(time_vec(CK_gt_0), C_D(CK_gt_0),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_D')
% subplot(3,1,3)
% semilogx(time_vec(CK_gt_0), C_K(CK_gt_0)./C_D(CK_gt_0),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_K/C_D')
% xlabel('t [s]')
% set(gcf,'color','w','position',[440   227   312   571])
% lh = legend('-dynamiclegend');
% set(lh,'location','northwest')
%
%
% figure(14)
% subplot(1,2,1)
% yyaxis left
% semilogx(time_vec(1:end-1),cumsum(C_K.*diff(time_vec)));hold on
% ylabel(' C_K')
% yyaxis right
% semilogx(time_vec(1:end-1),cumsum(C_D.*diff(time_vec)));hold on
% ylabel(' C_D')
% xlabel('t [s]')
%
% figure(4)
% subplot(1,2,1)
% yyaxis left
% semilogx(time_vec(1:end-1),cumsum(C_K.*diff(time_vec)./time_vec(2:end)));hold on
% ylabel(' C_K')
% yyaxis right
% semilogx(time_vec(1:end-1),cumsum(C_D.*diff(time_vec)./time_vec(2:end)));hold on
% ylabel(' C_D')
% xlabel('t [s]')
% subplot(1,2,2)
% semilogx(time_vec(1:end-1),cumsum(C_K.*diff(time_vec)./time_vec(2:end))./cumsum(C_D.*diff(time_vec)./time_vec(2:end)),'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% semilogx(time_vec(1:end-1),cumsum(C_K.*diff(time_vec)./time_vec(2:end))./C_D,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% legend('-dynamiclegend')
% CK_gt_0 = C_K_int_num>0;
% CK_gt_0 = true(length(C_K_int_num),1);

% figure(15)
% subplot(3,1,1)
% semilogx(time_vec(CK_gt_0), C_K_int_num(CK_gt_0)./(-C_K_int_denom(CK_gt_0)),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_K')
% set(gca,'xlim',[1e-3 30])
% set(gca,'fontsize',20)
% subplot(3,1,2)
% semilogx(time_vec(CK_gt_0), C_D(CK_gt_0),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_D')
% set(gca,'xlim',[1e-3 30])
% set(gca,'fontsize',20)
% subplot(3,1,3)
% semilogx(time_vec(CK_gt_0), (C_K_int_num(CK_gt_0)./(-C_K_int_denom(CK_gt_0)))./C_D(CK_gt_0),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_K/C_D')
% xlabel('t [s]')
% set(gca,'xlim',[1e-3 30])
% set(gcf,'color','w','position',[ 265     1   596   804])
% % legend('-dynamiclegend')
% set(gca,'fontsize',20)
% set(lh,'location','eastoutside','fontsize',12)

% C_K_int_denom_int = cumsum(C_K_int_denom.*diff(time_vec)./time_vec(2:end));
% figure(16)
% subplot(3,1,1)
% semilogx(time_vec(CK_gt_0), C_K_int_num(CK_gt_0)./(-C_K_int_denom_int(CK_gt_0)),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_K')
% subplot(3,1,2)
% semilogx(time_vec(CK_gt_0), C_D(CK_gt_0),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_D')
% subplot(3,1,3)
% semilogx(time_vec(CK_gt_0), (C_K_int_num(CK_gt_0)./(-C_K_int_denom_int(CK_gt_0)))./C_D(CK_gt_0),'linewidth',3,'displayname',sprintf('r_0=% 3d %s m',round(ic.r_0*1e6),char(956)));hold on
% title('C_K/C_D')
% xlabel('t [s]')
% set(gcf,'color','w')
% legend('-dynamiclegend')

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















