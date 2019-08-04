% prograss function

function [] = progress_fxn(status)
if status > 0.25 && status < 0.255
        fprintf('\n25 %% of integration complete\n')
elseif status > 0.5 && status < 0.505
        fprintf('50 %% of integration complete\n')
elseif status > 0.75 && status < 0.755
        fprintf('75 %% of integration complete\n')
end


