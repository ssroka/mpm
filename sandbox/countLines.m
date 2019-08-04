totLines = 0;
list_dir = {'inputFiles','mm_integration','mm_vars','plotScripts','util','test_MPM','approximationEndpoints'};
for d = 1:length(list_dir)
    currentDir = list_dir{d};
    currentFiles = dir(currentDir);
    files_list = {currentFiles.name};
    for f = 3:length(files_list)
        if strcmp(files_list{f}(end-1:end),'.m')
            file_name = strcat(list_dir{d},'/',files_list{f});
            fprintf('%s :',files_list{f})
            fid = fopen(file_name);
            numLine = 0;
            while ~feof(fid)
                junk_line =fgetl(fid);
                numLine = numLine + 1;
            end
            fprintf(' %d\n',numLine)
            totLines = totLines + numLine;
            fclose(fid);
        end
    end
end

list_indiv_files = {'integrate_Qs_QL_driver.m',...
                    'integrate_Qs_QL.m',...
                    'dropletEvolution.m',...
                    'dropletEvolution_driver.m',...
                    'rm_ASF_paths.m'};
                    

    for f = 3:length(list_indiv_files)
        if strcmp(list_indiv_files{f}(end-1:end),'.m')
            file_name = list_indiv_files{f};
            fprintf('%s :',list_indiv_files{f})
            fid = fopen(file_name);
            numLine = 0;
            while ~feof(fid)
                junk_line =fgetl(fid);
                numLine = numLine + 1;
            end
            fprintf(' %d\n',numLine)
            totLines = totLines + numLine;
            fclose(fid);
        end
    end



fprintf('total lines: %d\n\n',totLines)