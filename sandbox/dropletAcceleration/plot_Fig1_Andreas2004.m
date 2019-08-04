

clear;close all;clc

d = dir(sprintf('%s/..',pwd));
for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        addpath(sprintf('%s/../%s',pwd,d(ii).name))
    end
end

ic.p0 = 100000;
ic.T_a = 18;
ic.T_s_0 = 20;
ic.S0 = 34;
U10 = 10:10:50;%
r_0 = logspace(0,log10(4000))*1e-6;%
% Newton-Raphson for r_eq
maxIt = 1000; % the maxium number of iterations before aborting the iterative scheme
maxEr_uf  = calcVterm_sphere(max(r_0))*0.0001;  % the max error in u_f

tau_ac = zeros(length(r_0),length(U10));
tau_f = zeros(length(r_0),length(U10));
r_0_vec = zeros(length(r_0),length(U10));
C_D = zeros(length(r_0),length(U10));
tau_tot = zeros(length(r_0),length(U10));

for jj = 1:length(U10)
    ic.U10 = U10(jj);
    for ii = 1:length(r_0)
        ic.r_0 = r_0(ii);
        ic.m_s = compute_initial_drop_conditions(ic.T_s_0,ic.S0,ic.r_0,ic.p0);
        
        helpingAnonFxns;
        if ic.r_0 < 100*1e-6
            t_max = 0.5; % seconds
        elseif ic.r_0 < 500*1e-6
            t_max = 3; % seconds
        elseif ic.r_0 < 1000*1e-6
            t_max = 10; % seconds
        else
            t_max = 80; % seconds
        end
        
        u0    = 0; % m/s starting from rest
        
        ra = rho_a(Celsius2Kelvin(ic.T_a),ic.p0);
        rs = rho_s(ic.T_s_0,ic.r_0,ic.m_s,ic.S0/1000,ic.p0);
        options = odeset('AbsTol',1e-9,'RelTol',1e-9);
        [acel_time_vec, u_vec] = ode45(@(t,u) compute_dudt(t,u,ra,rs,ic.U10,ic.r_0,ic.T_a),[0 t_max],u0,options);
        
        tau_ac(ii,jj) = acel_time_vec(max(find(u_vec<ic.U10-ic.U10*0.623)));
        tau_tot(ii,jj) = acel_time_vec(max(find(u_vec<ic.U10-ic.U10*0.01)));
        r_0_vec(ii,jj) = ic.r_0;
        tau_f(ii,jj) = compute_tauf(ic.U10,ic.T_s_0,ic.r_0,ic.m_s,ic.S0/1000,ic.p0,ic.T_a,maxEr_uf,maxIt);
        C_D(ii,jj) = (rs*4/3*ic.r_0*trapz(acel_time_vec(1:end-1),diff(u_vec)./diff(acel_time_vec)))/(ra*ic.U10^2*min(tau_tot(ii,jj),tau_f(ii,jj)));
    end
end


%%
figure
LS = {'--',':','-',':','-.'};
lw = [1,1,1,2,1];
for kk = 1:length(U10)
yyaxis left
loglog(r_0_vec(:,kk)*1e6,tau_tot(:,kk),'linestyle',LS{kk},'color','g','marker','none','linewidth',lw(kk))
hold on
loglog(r_0_vec(:,kk)*1e6,tau_ac(:,kk),'linestyle',LS{kk},'color','b','marker','none','linewidth',lw(kk))
h1(kk) = loglog(r_0_vec(:,kk)*1e6,tau_f(:,kk),'linestyle',LS{kk},'color','k','marker','none','linewidth',lw(kk));
yyaxis right
loglog(r_0_vec(:,kk)*1e6,C_D(:,kk),'linestyle',LS{kk},'color','r','marker','none','linewidth',lw(kk))
end
legend(h1,'U_{10} = 10','U_{10} = 20','U_{10} = 30','U_{10} = 40','U_{10} = 50','location','west')
set(gcf,'color','w')
xlabel('r_0 [\mu m]')
yyaxis left
ylabel('\tau_{f},\tau_{tot},\tau_{ac} [s]')
set(gca,'xlim',[0.1 10000],'ylim',[1e-6 1e6],'ycolor','k')
t1 = text(10,40,'\tau_f');
t2 = text(10,4e-2,'\tau_{tot}');
t3 = text(10,9e-5,'\tau_{ac}');
t4 = text(0.15,200,'T_a=T_s=20^{\circ}C');
set([t1 t2 t3 t4],'fontsize',16)
set(t2,'color','g')
set(t3,'color','b')

yyaxis right
ylabel('C_D')
set(gca,'ycolor','r')
cdtxt = text(200,0.5,'C_D');
set(cdtxt,'color','r','fontsize',16)
title('Figure 1 Andreas 2004')
set(gcf,'position',[ 86   421   484   368])
%{

addpath ~/Documents/MATLAB/util/export_fig/
export_fig(gcf,'tau_CD.pdf','-m2')
copyfile('tau_CD.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

    





%}



%%

run /Users/ssroka/Documents/MATLAB/ASF/calcCoeff_1drop/sandbox/rm_ASF_paths












