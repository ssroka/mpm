
clear
close all
clc

% print the RH's, DT's, and Ts's
% RH_list = [88.0, 90.0, 92.0, 94.0, 96.0, 98.0];
RH_list = [80.0, 82.0, 84.0, 86.0]
DT_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
Ts_list = [27.0, 28.0, 29.0]

for i_RH = 1:2;%(length(RH_list))	
 for i_DT = 1:(length(DT_list))	
  for i_Ts = 1:(length(Ts_list))	
   % exp_num = 1000*(1+i_Ts) + 100*(i_DT-1) + 10*(i_RH)
   % exp_num = 100*(1+i_Ts) + 10*(i_DT-1) + 1*(i_RH)
   exp_num = 1000*(4+i_Ts) + 100*(i_DT-1) + 10*(i_RH)
   file_str = num2str(exp_num);
   dig3 = str2double(file_str(3));
   if dig3<1 || dig3>6
      continue
   end
   res_dir = sprintf('~/mpm/sandbox/results/Exp%d',exp_num);
   store_rt_timehistories(res_dir);
  end
 end
end


 
