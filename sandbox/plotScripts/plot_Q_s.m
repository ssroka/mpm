% plot Q_s
% clear; clc;
%resultsDir = {'/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/Andreas08_Fig2_in_results',...
%'/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/Andreas08_Fig2_in_results_testLargeIntError'};
% resultsDir = {'/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/Andreas08_Fig2_in_results_test1NayayFxns'};
function [] = plot_Q_s(resultsDir)

n = 0;
for jj = 1:length(resultsDir)
    if ~strcmp(resultsDir{jj}(1:6),'/Users')
        resultsDir{jj} = strcat(pwd,'/',resultsDir{jj});
    end
mydir = dir(resultsDir{jj});
str = {mydir.name};
pat = 'r_0_\w*\.mat';
files_tmp = regexp(str,pat,'match');
for ii = 1:length(files_tmp)
  if ~isempty(files_tmp{ii})
    files{ii+n}{1} = strcat(resultsDir{jj},'/',files_tmp{ii}{1});
  end
end
  n = n + length(files_tmp);
end
% thisFile = mfilename('fullpath');
% cd(thisFile(1:max(strfind(thisFile,'/'))))

% --------------------------------------------------------
%                          Qs
% --------------------------------------------------------

count = 1;
for ii = 1:length(files)
    if ~isempty(files{ii})
        load(files{ii}{1})
        i_star = sum(time_vec<tau_f);
        if i_star==length(time_vec)
            T_at_tau_f = T_eq_fmm;
            r_at_tau_f = r_eq_fmm;
        else
            T_at_tau_f = (T_s_t(i_star+1)-T_s_t(i_star))/(time_vec(i_star+1)-time_vec(i_star))*(tau_f - time_vec(i_star))+T_s_t(i_star);
            r_at_tau_f = (r_t(i_star+1)-r_t(i_star))/(time_vec(i_star+1)-time_vec(i_star))*(tau_f - time_vec(i_star))+r_t(i_star);
        end
        if ic.r_0>85*1e-6
         dFdr0  = SSGF_OrtizSuslow(ic.r_0*1e6);% argument comes in as micrometers
    New_Vol_flux = 4/3*pi*ic.r_0^3*dFdr0; % units are m^2/(s * micrometers)

        
%         dFdr0  = SSGF_Fairall94(r_t(1)*1e6);% argument comes in as micrometers
%         W_u = 3.8*1e-6*ic.U10^(3.4); % whitecap fraction
%         New_Vol_flux = dFdr0*W_u;
        
        
        s_at_tau_f = ic.m_s./(ic.m_s+compute_m_w(T_at_tau_f,r_at_tau_f,ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt));
        Q_L(count) = rho_s(T_s_t(1),r_t(1),ic.m_s,ic.S0/1000,ic.p0).*L_v(T_s_t(1),ic.S0).*(1-(r_at_tau_f/r_t(1))^3).*New_Vol_flux;
        Q_L_tau_f(count) = rho_s(T_s_t(1),r_t(1),ic.m_s,ic.S0/1000,ic.p0).*L_v(T_at_tau_f,ic.S0).*(1-(r_at_tau_f/r_t(1))^3).*New_Vol_flux;
        Q_L_eq(count) = rho_s(T_s_t(1),r_t(1),ic.m_s,ic.S0/1000,ic.p0).*L_v(T_s_t(1),ic.S0).*(1-(r_eq_fmm/r_t(1))^3).*New_Vol_flux;
        Q_s(count) = rho_s(T_s_t(1),r_t(1),ic.m_s,ic.S0/1000,ic.p0).*c_ps(T_s_t(1),ic.S0,ic.p0).*(T_s_t(1)-T_eq_fmm).*(1-exp(-tau_f/tau_T_fmm)).*New_Vol_flux;
        Q_s_r_0(count) = r_t(1);
        tau_f_vec(count) = tau_f;
        VolFluxes(count) = dFdr0;
        tau_T_vec(count) = tau_T_fmm;
        rho_s_vec(count) = rho_s(T_s_t(1),r_t(1),ic.m_s,ic.S0/1000,ic.p0);
        count = count + 1;
        end
    end
end

h = loglog(Q_s_r_0*1e6,Q_s','*','color',[0 0 1],'displayname','Q_s');
hold on
loglog(Q_s_r_0*1e6,Q_L','r*','displayname','Q_L')
% loglog(Q_s_r_0*1e6,Q_L_tau_f','kd')
% Q_s_fromWebPlotDig;
% loglog(Q_s_fromWPD(:,1),Q_s_fromWPD(:,2),'bo')
% loglog(Q_L_fromWPD(:,1),Q_L_fromWPD(:,2),'ko')
% set(gca,'xlim',[0.1 1000],'ylim',[1e-7 1e2])

set(gca,'xlim',[1 1e4],'ylim',[1e-6 1e0])
title('Q_s and Q_L')
legend('-dynamiclegend')
set(gcf,'color','w')
xlabel('r_0 [\mu m]')
ylabel('$$Q_s, Q_L \qquad [W m^{-2} \mu m^{-1}]$$','interpreter','latex')
mystr = sprintf('SSGF from Ortiz-Suslow 2016 Fig 8');
text(2,.10,mystr)
SSGFText = text(2,.056,'$$U_{10} = 36 m/s$$');
set(SSGFText,'interpreter','latex')
% loglog([99.33862756760169,99.33862756760169,99.33862756760169,99.33862756760169],[0.03012947724726992,0.027885481717262916,0.024828928805603275,0.022979709690469627],'r*')
% loglog([100.6657756892686, 100.6657756892686, 100.6657756892686, 100.6657756892686, 100.6657756892686],[0.15910384392721516,0.17190722018585747,0.185740907463857,0.20068781676649666,0.2168375310987433],'r*')


figure
loglog(Q_s_r_0*1e6,tau_f_vec./(0.015*ic.U10^2),'*');
title('tau f from Andreas 1990')
xlabel('r_0')
ylabel('tau_f')

return 

%{
    addpath ~/Documents/MATLAB/util/export_fig/
    
 export_fig(gcf,'Q_s_Q_L_Suslow.pdf','-m2')
 copyfile('Q_s_Q_L_Suslow.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

    
    %}

figure
loglog(Q_s_r_0*1e6,tau_f_vec,'*');
title('tau f')
xlabel('r_0')
ylabel('tau_f')

figure
loglog(Q_s_r_0*1e6,tau_T_vec,'*');
title('tau T')
xlabel('r_0')
ylabel('tau_T')


figure
loglog(Q_s_r_0*1e6,VolFluxes,'k*');
FairallWebPlotDig_Fig3;
hold on 
loglog(Fairall_Fig3_WebPlotDig_data(:,1),1e-14*Fairall_Fig3_WebPlotDig_data(:,2),'bo');
title('VolFluxes ')
xlabel('r_0')
ylabel('(4 \pi r_0^3 /3 ) dF/dr_0')

figure
loglog(Q_s_r_0*1e6,0.015*ic.U10^2./tau_f_vec,'*');
title('u f ')
xlabel('r_0')
ylabel('u_f')

end







