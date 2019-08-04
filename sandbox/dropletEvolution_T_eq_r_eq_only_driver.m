function [] = dropletEvolution_T_eq_r_eq_only_driver(slurmJobID)
%% User Input

% define initial microphycial constants and parameters
inputFile = 'Exp20_in'; % Andreas05_Fig11_in Andreas08_Fig2_in Andreas90_Fig09_in CommitteeMeeting_500micormeterdrop


%% ----------------------------- Set Parameters from runfile  --------------------------
filetext = fileread(sprintf('inputFiles/%s.m',inputFile));

%  Set Nayar Function Flag 
global Nayar_flag ;
start_char = strfind(filetext,'Nayar_flag');
last_char = min(strfind(filetext(start_char:start_char+100),sprintf('\n')));
eval(filetext(start_char:start_char+last_char))

% Set Save Script 
start_char = strfind(filetext,'saveScript');
last_char = min(strfind(filetext(start_char:end),sprintf('\n')));
eval(filetext(start_char:start_char+last_char))

% Set Save Directory 
start_char = strfind(filetext,'saveDir');
last_char = min(strfind(filetext(start_char:start_char+100),sprintf('\n')));
eval(filetext(start_char:start_char+last_char))
%% ----------------------------- Create Entry in Log File --------------------------
fid = fopen('logfile.txt','a');
fprintf(fid,'-------------------------\n');
if nargin==1; fprintf(fid,'Job ID: % 9.0f\n',slurmJobID); end
fprintf(fid,'Run Name: %s  \n',saveDir);
logStr = {'off','on'};
fprintf('Nayar Functions: %s\n',logStr{Nayar_flag + 1});
fprintf(fid,'Nayar Functions: %s\n',logStr{Nayar_flag + 1});
fclose(fid);

%% ----------------------------- setup paths --------------------------
pwd_tmp = pwd;
dirstr = pwd_tmp(1:strfind(pwd,'sandbox')+length('sandbox')-1);
d = dir(dirstr);
addpath(dirstr)
for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        addpath(sprintf('%s/%s',dirstr,d(ii).name))
    end
end

%% ----------------------------- begin droplet evolution --------------------------
dropletEvolution_T_eq_r_eq_only;


%% ----------------------------- email results --------------------------
% setpref('Internet','SMTP_Server','imap.exchange.mit.edu')
% sendmail('ssroka@mit.edu','Job finished running','Your results are here',...
%          sprintf('/Users/ssroka/Documents/MATLAB/ASF/calcCoeff_1drop/sandbox/results/%s/CK_CD_rat_plot.png',saveDir))

%% ----------------------------- remove user-added paths --------------------------
rm_ASF_paths;
