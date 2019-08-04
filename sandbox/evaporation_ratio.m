clear;close all;clc
[CK1, CD1] = calc_CK_CD_1_drop(strcat(pwd,'/results/Exp1_01/r_0_100_Exp1_01.mat'));

rat(1,1) = CK1/CD1;

% [CK2, CD2] = calc_CK_CD_1_drop(strcat(pwd,'/results/Exp2/r_0_100_Exp2.mat'));
% 
rat(1,2) = CK1/CD1;

rat_exact =  linspace(0.5,1)';

percentEvap =repmat(rat_exact,1,2)*100./repmat(rat,length(rat_exact),1);


plot(rat_exact,percentEvap(:,1),'b--','linewidth',2,'displayname',sprintf('U_{10} = 36, C_K/C_D = %1.2f',rat(1)))
hold on
plot(rat_exact,percentEvap(:,2),'r-','linewidth',2,'displayname',sprintf('U_{10} = 20, C_K/C_D = %1.2f',rat(2)))


xl = get(gca,'xlim');
set(gca,'ylim',[0 100])
yl = get(gca,'ylim');

text((xl(2)-xl(1))*.2+xl(1),(yl(2)-yl(1))*.8+yl(1),'evaporate')
text((xl(2)-xl(1))*.7+xl(1),(yl(2)-yl(1))*1/6+yl(1),'re-enter')
ylabel('100 drops')
xlabel('desired C_K/C_D')
title('Evaporation Percentages: r_0 = 100 \mu m')
set(gcf,'color','w')

lh = legend('-dynamiclegend');
set(lh,'location','eastoutside')

% c_pm = 1017;

%{
    addpath ~/Documents/MATLAB/util/export_fig/
 export_fig(gcf,'Evaporation_Percent.pdf','-m2')
 copyfile('Evaporation_Percent.pdf','~/MIT/Research/EmanuelGroup/PhD_ThesisProposal_abstract/AgencyProposal/')

    
%}





