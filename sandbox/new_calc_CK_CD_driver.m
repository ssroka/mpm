clear;close all;clc
save_new_data_flag = false;

% file_pre = '/Users/ssroka/Documents/MATLAB/ASF/calcCoeff_1drop/sandbox/results/';
file_pre = '/Users/ssroka/Documents/MATLAB/ASF/calcCoeff_1drop/fromEngaging/Exp1_100_101_104/';
delete fig21.mat
for i = 2:2:12
    if i<=2
        file_suf = sprintf('r_0_%d_Exp1_104.mat',25*i);
    elseif i<=5
        file_suf = sprintf('r_0_%d_Exp1_101.mat',25*i);
    else
        file_suf = sprintf('r_0_%d_Exp1_100.mat',25*i);
    end


file_name = [file_pre file_suf];

%% load variables

load(file_name)

% yyaxis left
% semilogx(time_vec,T_s_t,'*')
% ylabel('T [deg C]')
% hold on
% yyaxis right
% semilogx(time_vec,r_t,'*')
% ylabel('r [m]')
% xlabel('t [s]')

%% integrate Qs and QL

new_calc_CK_CD(file_name,save_new_data_flag);

end
