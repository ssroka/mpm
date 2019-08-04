%clear;clc;close all
function [] = driverFromSlurm(slurmJobID)
%% User Input

% define initial microphycial constants and parameters
inputFile = 'Exp1_in'; % Andreas05_Fig11_in Andreas08_Fig2_in Andreas90_Fig09_in CommitteeMeeting_500micormeterdrop


%% ----------------------------- Nayar Functions --------------------------
global Nayar_flag ;


%% ----------------------------- Wrtie to Log File  --------------------------

if nargin>0
fid =	fopen('logFile.txt','at');
fprintf(fid,'\n----------------------------\n');
fprintf(fid,'Job ID: % 9f\n',slurmJobID);
fprintf(fid,'input file : %s\n',inputFile);
logStr = {'false','true'};
fprintf(fid,'Nayar flag: %s\n',logStr{Nayar_flag+1});
fclose(fid);
end


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
%% ----------------------------- remove user-added paths --------------------------
d = dir(pwd);

for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        rmpath(sprintf('%s/%s',pwd,d(ii).name))
    end
end
