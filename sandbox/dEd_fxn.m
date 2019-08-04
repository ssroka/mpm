function [] = dropletEvolution_driver(slurmJobID)
%% User Input

% define initial microphycial constants and parameters
inputFile = 'Exp1_in'; % Andreas05_Fig11_in Andreas08_Fig2_in Andreas90_Fig09_in CommitteeMeeting_500micormeterdrop


%% ----------------------------- Set Parameters from runfile  --------------------------
filetext = fileread(sprintf('inputFiles/%s.m',inputFile));

%  Set Nayar Function Flag 
global Nayar_flag ;
start_char = strfind(filetext,'Nayar_flag');
last_char = min(strfind(filetext(start_char:start_char+100),sprintf('\n')));
eval(filetext(start_char:start_char+last_char))

% Set Save Directory 
start_char = strfind(filetext,'saveDir');
last_char = min(strfind(filetext(start_char:start_char+100),sprintf('\n')));
eval(filetext(start_char:start_char+last_char))

%% ----------------------------- Create Entry in Log File --------------------------
fid = fopen('logfile.txt','a');
fprintf(fid,'-------------------------\n');
fprintf(fid,'Job ID: % 9.0f\n',slurmJobID);
fprintf(fid,'Run Name: %s  \n',saveDir);
logStr = {'off','on'};
fprintf('Nayar Functions: \n',logStr{Nayar_flag + 1});
fprintf(fid,'Nayar Functions: \n',logStr{Nayar_flag + 1});
fclose(fid);

%% ----------------------------- setup paths --------------------------

d = dir(pwd);
for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        addpath(sprintf('%s/%s',pwd,d(ii).name))
    end
end

%% ----------------------------- begin droplet evolution --------------------------
dropletEvolution;

%% ----------------------------- compute CK-CD --------------------------
plot_CK_CD_ratio({sprintf('results/%s',saveDir)});
savefig(gcf,sprintf('results/%s/CK_CD_rat_plot.fig',saveDir))
print(gcf,sprintf('results/%s/CK_CD_rat_plot',saveDir),'-dpng')

%% ----------------------------- email results --------------------------
setpref('Internet','E_mail','ssroka@mit.edu')
sendmail('ssroka@mit.edu','Job finished running','Your results are here',...
         sprintf('/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/results/%s/CK_CD_rat_plot.png',saveDir))

%% ----------------------------- remove user-added paths --------------------------
d = dir(pwd);

for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        rmpath(sprintf('%s/%s',pwd,d(ii).name))
    end
end
