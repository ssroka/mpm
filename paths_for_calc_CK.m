paths_2_add = {
    '~/mpm/sandbox/util';
    '~/util';
    '~/mpm/sandbox/mm_vars/';
    '~/mpm/sandbox/Q_s/';
    '~/mpm/sandbox/NayarFxns/';
    '~/util/distinguishable_colors';
    '~/util/othercolor';
    '~/mpm_figures/';
    '~/mpm/SGFs';
    '~/mpm/sandbox/dropletAcceleration'
    };
switch action_str
    case 'add'
        for ipath = 1:length(paths_2_add)
            addpath(paths_2_add{ipath})
        end
    case 'remove'
        for ipath = 1:length(paths_2_add)
            rmpath(paths_2_add{ipath})
        end
    otherwise
        error('Please select a valid path action, \n''add'' or ''remove''')
end


