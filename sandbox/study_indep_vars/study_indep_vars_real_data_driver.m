





clear;close all;clc
save_new_data_flag = true;

file_pre = '/Users/ssroka/Documents/MATLAB/ASF/calcCoeff_1drop/sandbox/study_indep_vars/';
delete fig21.mat
for i = 2
    for delta_T = 1:2
        switch delta_T
            case 1
                file_suf = sprintf('r_0_%d_Exp5_100_3.mat',50*i);
            case 2
                file_suf = sprintf('r_0_%d_Exp5_100_3.mat',50*i);
            case 3
                file_suf = sprintf('r_0_%d_Exp5_100_3.mat',50*i);
            case 4
                file_suf = sprintf('r_0_%d_Exp5_100_3.mat',50*i);
            case 5
                file_suf = sprintf('r_0_%d_Exp5_100_3.mat',50*i);
        end
        
        file_name = [file_pre file_suf];
        
        %% load variables
        load(file_name)
        U10s = [10:10:50];
        for j = 1:length(U10s)
            if j == 1
                study_indep_vars_real_data(file_name,save_new_data_flag,U10s(j),true,delta_T);
            else
                study_indep_vars_real_data(file_name,save_new_data_flag,U10s(j),true,delta_T);
            end
        end
    end
end






