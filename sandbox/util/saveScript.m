% all results go here
if ~exist('results','dir')
    mkdir('results')
end
% results from this specific simulation go here
if ~exist('saveDir','var')
    saveDir = sprintf('%s_%s',inputFile,datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
end
if ~exist(sprintf('results/%s/',saveDir),'dir')
    mkdir(sprintf('results/%s/',saveDir))
end

% results from this iteration of this simulation go here
results_loc = sprintf('results/%s/r_0_%s_%s',saveDir,strrep(num2str(round(r_0*1e6,4)),'.','_'),saveDir);

% save these results
save(results_loc,'ic','time_vec','r_t','T_s_t','T_eq','T_eq_fmm','r_eq','r_eq_fmm','tau_T','tau_T_fmm','tau_r','tau_r_fmm','tau_f','t_final')





















