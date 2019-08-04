

% resultsDir = {'/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/results/Exp1_01'};
resultsDir = {'/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/fromEngaging/Exp1_2'};

[r_0_vec_tmp,CK_CD_ratio_tmp,CK_vec_tmp,CD_vec_tmp] = plot_CK_CD_ratio(resultsDir);

[r_0_vec, ind] = sort(r_0_vec_tmp);
CK_CD_ratio = CK_CD_ratio_tmp(ind);
CK_vec = CK_vec_tmp(ind);
CD_vec = CD_vec_tmp(ind);

fo = fitoptions('Method','NonlinearLeastSquares');
f1 = fit(r_0_vec*1e6,CK_CD_ratio,'exp1',fo);
f2 = fit(r_0_vec*1e6,CK_vec,'exp1',fo);
f3 = fit(r_0_vec*1e6,CD_vec,'exp1',fo);
%% 
close all
figure
% subplot(4,6,[3:6,9:12,15:18,21:24])
subplot(3,1,3)
semilogx(r_0_vec*1e6,CK_CD_ratio,'*')
hold on
f = fit(r_0_vec*1e6,CK_CD_ratio,'exp1');
semilogx(r_0_vec*1e6,f.a*exp(r_0_vec*1e6*f.b),'r-','linewidth',1)
% text(60,120,sprintf('Fit: %3.1fe^{%1.3fr_0}',f.a,f.b))
xlabel('r_0 [\mu m]')
title(strcat('C_K/C_D','    ',sprintf('    Fit: %3.1fe^{%1.3fr_0}',f.a,f.b)))

% subplot(4,6,[1:2,7:8])
subplot(3,1,1)
semilogx(r_0_vec*1e6,CK_vec,'*')
hold on
f = fit(r_0_vec*1e6,CK_vec,'exp1');
semilogx(r_0_vec*1e6,f.a*exp(r_0_vec*1e6*f.b),'r-','linewidth',1)
% text(60,3,sprintf('Fit: %3.1fe^{%1.3fr_0}',f.a,f.b))
title(strcat('C_K','    ',sprintf('    Fit: %3.1fe^{%1.3fr_0}',f.a,f.b)))


% subplot(4,6,[13:14,19:20])
subplot(3,1,2)
semilogx(r_0_vec*1e6,CD_vec,'*')
hold on
fe = fit(r_0_vec*1e6,CD_vec,'exp2');
semilogx(r_0_vec*1e6,fe.a*exp(r_0_vec*1e6*fe.b)+fe.c*exp(r_0_vec*1e6*fe.d),'r-','linewidth',1)
% text(60,0.025,sprintf('Fit: %1.2fe^{%1.3fr_0}+%1.2fe^{%1.3fr_0}',fe.a,fe.b,fe.c,fe.d))
title(strcat('C_D','    ',sprintf('    Fit: %1.2fe^{%1.3fr_0}+%1.2fe^{%1.3fr_0}',fe.a,fe.b,fe.c,fe.d)))
set(gcf,'color','w','position',[112   230   301   568])




%%
lrgdrops = 4:length(r_0_vec);
figure
% subplot(4,6,[3:6,9:12,15:18,21:24])
subplot(3,1,3)
semilogx(r_0_vec(lrgdrops)*1e6,CK_CD_ratio(lrgdrops),'*')
hold on
f = fit(r_0_vec(lrgdrops)*1e6,CK_CD_ratio(lrgdrops),'exp1');
semilogx(r_0_vec(lrgdrops)*1e6,f.a*exp(r_0_vec(lrgdrops)*1e6*f.b),'r-','linewidth',1)
% text(60,120,sprintf('Fit: %3.1fe^{%1.3fr_0}',f.a,f.b))
xlabel('r_0 [\mu m]')
title(strcat('C_K/C_D','    ',sprintf('    Fit: %3.1fe^{%1.3fr_0}',f.a,f.b)))

% subplot(4,6,[1:2,7:8])
subplot(3,1,1)
semilogx(r_0_vec(lrgdrops)*1e6,CK_vec(lrgdrops),'*')
hold on
f = fit(r_0_vec(lrgdrops)*1e6,CK_vec(lrgdrops),'exp1');
semilogx(r_0_vec(lrgdrops)*1e6,f.a*exp(r_0_vec(lrgdrops)*1e6*f.b),'r-','linewidth',1)
% text(60,3,sprintf('Fit: %3.1fe^{%1.3fr_0}',f.a,f.b))
title(strcat('C_K','    ',sprintf('    Fit: %3.1fe^{%1.3fr_0}',f.a,f.b)))


% subplot(4,6,[13:14,19:20])
subplot(3,1,2)
semilogx(r_0_vec(lrgdrops)*1e6,CD_vec(lrgdrops),'*')
hold on
fe = fit(r_0_vec(lrgdrops)*1e6,CD_vec(lrgdrops),'exp2');
semilogx(r_0_vec(lrgdrops)*1e6,fe.a*exp(r_0_vec(lrgdrops)*1e6*fe.b)+fe.c*exp(r_0_vec(lrgdrops)*1e6*fe.d),'r-','linewidth',1)
% text(60,0.025,sprintf('Fit: %1.2fe^{%1.3fr_0}+%1.2fe^{%1.3fr_0}',fe.a,fe.b,fe.c,fe.d))
title(strcat('C_D','    ',sprintf('    Fit: %1.2fe^{%1.3fr_0}+%1.2fe^{%1.3fr_0}',fe.a,fe.b,fe.c,fe.d)))
set(gcf,'color','w','position',[112   230   301   568])




%{

addpath ~/Documents/MATLAB/util/export_fig/
export_fig(gcf,'CD_CK_rat_zoom.pdf','-m2')
copyfile('CD_CK_rat_zoom.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

addpath ~/Documents/MATLAB/util/export_fig/
export_fig(gcf,'CD_CK_rat.pdf','-m2')
copyfile('CD_CK_rat.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

    





%}



















