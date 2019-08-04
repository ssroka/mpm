
clear
close all
clc



for file_num = 100:170
%	try
	file_str = num2str(file_num);
	dig3 = str2double(file_str(3));
	if dig3<1 || dig3>6
		continue
        end
res_dir = sprintf('~/mpm/sandbox/results/Exp%d',file_num);
store_rt_timehistories(res_dir);
%{ catch
	cd('~/mpm/sandbox/util')
	continue
%}end
end


 
