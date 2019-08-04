% plot CK and CD
function [] =plot_CK_CD(resultsFiles)
close all;clc
load(resultsFiles)
helpingAnonFxns;

[CK_tmp, C_D] = calc_CK_CD_1_drop(resultsFiles);
C_k = ones(length(C_D),1)*CK_tmp;
%% unscaled
figure
yyaxis left
semilogx(time_vec,C_k,'b--','displayname','C_K')
hold on
semilogx(time_vec,C_D,'b.-','displayname','C_D')
ylabel('C_K , C_D')
yyaxis right
semilogx(time_vec,C_k./C_D,'r','displayname','C_K/C_D')
ylabel('ratio C_K/C_D')
title('C_K , C_D, and C_K/C_D')
xlabel('t [s]')
legend2=legend('-dynamiclegend');
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
tau_f_line = line([tau_f tau_f],[ylim(1) ylim(2)]);
set(tau_f_line,'color','k')
set(gca,'ylim',ylim);

xtickmarks = get(gca,'xtick');
xticklabels_tmp = get(gca,'xticklabels');
[xtickmarks ind] = sort([xtickmarks(1:end-1),tau_f]);
xticklabels = [xticklabels_tmp(1:end-1); {'\tau_f'}];

text('FontSize',14,...
    'String',{'r_0 = 100 µm','T_s = 20 ºC','T_a = 18 ºC','RH = 90 %'},...
    'Position',[0.0167797628450504 1.25e4 0]);

% frac_up = 1/6;
% frac_right = 2.036/20;
% % see ASCII code table
% param_str= sprintf('r_0 = %g %sm \nT_s = %g %sC \nT_a = %g %sC\nRH = %g %s',ic.r_0*1e6,char(181),ic.T_s,char(186),ic.T_a,char(186),ic.RH,char(37));
% problemParam_str = text(log10((xlim(2)-xlim(1))*frac_right+(xlim(1))),(ylim(2)-ylim(1))*frac_up+ylim(1),param_str,'interpreter','tex');
% set(problemParam_str,'fontsize',14)
set(gca,'xtick',xtickmarks,'xticklabel',xticklabels(ind))
set(legend,'fontsize',14,'position',[ 0.14553      0.34691      0.15357      0.19524])
set(gcf,'color','w')
set(legend2,...
    'Position',[0.152672857142855 0.60167190476191 0.15357 0.195240000000002],...
    'FontSize',14);
%{
export_fig(gcf,'C_K_C_D.pdf','-m2')
copyfile('C_K_C_D.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

    %}


%% scaled
scale_factor = min(C_D)./max(C_k)
C_k_scale = scale_factor*C_k;

figure1 = figure('Color',[1 1 1]);

yyaxis left
semilogx(time_vec,C_k_scale,'b--','displayname','C_K^*')
hold on
semilogx(time_vec,C_D,'b.-','displayname','C_D')
ylabel('C_K^* , C_D')
yyaxis right
semilogx(time_vec,C_k_scale./C_D,'r','displayname','C_K^* /C_D')
ylabel('ratio C_K^* /C_D')
title('C_K^*, C_D, and C_K^* /C_D')
xlabel('t [s]')
legend1 = legend('-dynamiclegend');
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
tau_f_line = line([tau_f tau_f],[ylim(1) ylim(2)]);
set(tau_f_line,'color','k')
set(gca,'ylim',ylim);

xtickmarks = get(gca,'xtick');
xticklabels_tmp = get(gca,'xticklabels');
[xtickmarks int] = sort([xtickmarks,tau_f]);
xticklabels = [xticklabels_tmp; {'\tau_f'}];
xticklabels = xticklabels(int);
xticklabels(xtickmarks==10) = {'  '};

frac_up = 1/7;
frac_right = 2.036/20;
% see ASCII code table
% param_str= sprintf('r_0 = %g %sm \nT_s = %g %sC \nT_a = %g %sC\nRH = %g %s',ic.r_0*1e6,char(181),ic.T_s,char(186),ic.T_a,char(186),ic.RH,char(37));
% problemParam_str = text(log10((xlim(2)-xlim(1))*frac_right+(xlim(1))),(ylim(2)-ylim(1))*frac_up+ylim(1),param_str,'interpreter','tex');
% set(problemParam_str,'fontsize',14)
text('FontSize',14,...
    'String',{'r_0 = 100 µm','T_s = 20 ºC','T_a = 18 ºC','RH = 90 %'},...
    'Position',[0.0167797628450504 .5 0]);
set(gca,'xtick',xtickmarks,'xticklabel',xticklabels)
% set(legend,'fontsize',14,'position',[ 0.14553      0.33691      0.15357      0.19524])
set(gcf,'color','w')

set(legend1,...
    'Position',[0.152672857142855 0.60167190476191 0.15357 0.195240000000002],...
    'FontSize',14);
%{
    addpath ~/Documents/MATLAB/util/export_fig/
    
 export_fig(gcf,'C_K_C_D_scaled.pdf','-m2')
 copyfile('C_K_C_D_scaled.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

    
    %}




