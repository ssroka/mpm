helpingAnonFxns;
close all

%% full microphysical model
figure(1)
yyaxis left
semilogx(time_vec,T_s_t,'--b')
hold on
yyaxis right
semilogx(time_vec,1e6*r_t,'--r')

%% estimates microphysical model
% estimates have different support in time than in full model
time_vec_for_T_est = time_vec(time_vec<=350);
time_vec_for_r_est = [time_vec;logspace(log10(time_vec(end)),log10(4000))'];

yyaxis left
T_est = Kelvin2Celsius(T_eq) + (ic.T_s - Kelvin2Celsius(T_eq))*exp(-time_vec_for_T_est/tau_T);
semilogx(time_vec_for_T_est,T_est ,'-b')
ylabel('Droplet Temperature [^{\circ} C]')
ylim([16 21])
yyaxis right
r_est = r_eq + (ic.r_0 - r_eq)*exp(-time_vec_for_r_est/tau_r);
semilogx(time_vec_for_r_est,1e6*r_est,'-r')
ylabel('Droplet Radius (\mu m)')
ylim([50 110])
xlim([0.01 5000])
xlabel('Time Since Formation (s)')
params_on_graph ={' r_0 = 100 \mu m ',' T_s = 20^{\circ}C ',' T_a = 18^{\circ}C ',' RH = 90%'};
text(200,100,params_on_graph)
title('Updated Meterological Functions')
set(gcf,'color','w')
yyaxis left
text(0.06,16.75,'Full Model')
text(0.06,16.5,'Exponential Approximation')
legend_Line1 = line([0.015 0.05],16.75*[1 1]);
set(legend_Line1,'color','k','linewidth',1)
legend_Line2 = line([0.015 0.05],16.5*[1 1]);
set(legend_Line2,'color','k','linewidth',1,'linestyle','--')

%{
    addpath ~/Documents/MATLAB/util/export_fig/
    
 export_fig(gcf,'r_T_updated.pdf','-m2')
 copyfile('r_T_updated.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

    
    %}
