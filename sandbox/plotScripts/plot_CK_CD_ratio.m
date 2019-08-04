function [r_0_vec,CK_CD_ratio,CK_vec,CD_vec] = plot_CK_CD_ratio(resultsDir)


n = 0;
for jj = 1:length(resultsDir)
    if ~strcmp(resultsDir{jj}(1:6),'/Users')
        resultsDir{jj} = strcat(pwd,'/',resultsDir{jj});
    end
    mydir = dir(resultsDir{jj});
    str = {mydir.name};
    pat = 'r_0_\w*\.mat';
    files_tmp = regexp(str,pat,'match');
    for ii = 1:length(files_tmp)
        if ~isempty(files_tmp{ii})
            files{ii+n}{1} = strcat(resultsDir{jj},'/',files_tmp{ii}{1});
        end
    end
    n = n + length(files_tmp);
end
added_path_flag = 0;
if isempty(which('compute_Teq.m'))
addpath /Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/
added_path_flag = 1;
end
nonEmptyCells = find(~cellfun(@isempty,files));
CK_CD_ratio = zeros(length(nonEmptyCells),1);
r_0_vec = zeros(length(nonEmptyCells),1);
CK_vec = zeros(length(nonEmptyCells),1);
CD_vec = zeros(length(nonEmptyCells),1);
for ii = 1:length(nonEmptyCells)
        [CK_vec(ii), CD_vec(ii)] = calc_CK_CD_1_drop(files{nonEmptyCells(ii)}{1});
        CK_CD_ratio(ii) = CK_vec(ii)/CD_vec(ii);
        load(files{nonEmptyCells(ii)}{1},'ic')
        r_0_vec(ii) = ic.r_0;
        fprintf('r_0 = % 3.1f\n',ic.r_0*1e6)
end



semilogx(r_0_vec*1e6,CK_CD_ratio,'*')
xlabel('r_0 [\mu m]')
ylabel('CK/CD')
figure
plot(r_0_vec*1e6,CK_CD_ratio,'*')
xlabel('r_0 [\mu m]')
ylabel('CK/CD')

if added_path_flag
rmpath /Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/
end





