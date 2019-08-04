clear;close all;clc


file_pre = '/Users/ssroka/Documents/MATLAB/ASF/calcCoeff_1drop/sandbox/results/';
for i = 1:4
file_suf = sprintf('Exp1_4/r_0_%d_Exp1_4.mat',50*i);

file_name = [file_pre file_suf];

%% load variables

load(file_name)

yyaxis left
semilogx(time_vec,T_s_t,'*')
ylabel('T [deg C]')
hold on
yyaxis right
semilogx(time_vec,r_t,'*')
ylabel('r [m]')

xlabel('t [s]')

%% integrate Qs and QL

[Qs,QL,t] = integrate_Qs_QL(file_name);

end











