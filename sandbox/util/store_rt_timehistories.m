function [] = store_rt_timehistories(directory_to_save,filename)
my_dir = pwd;
copyfile('mpm_sim.m',directory_to_save) 
cd(directory_to_save)
d = dir;
d_n = {d.name};
nc = 1;
for id = 1:length(d_n)
	if length(d_n{id})>4 && strcmp(d_n{id}(end-2:end),'mat')
		d_mat{nc} = d_n{id};
		nc = nc + 1;
	end
end
count = 1;
for im = 1:length(d_mat)
      fprintf('%s\n',d_mat{im});
      load(sprintf('%s',d_mat{im}));
      whos
      ic.r_0
      t_final
      timehistory(count) = mpm_sim(d_mat{im},r_t,T_s_t,time_vec,ic);
      count = count + 1;
end
if nargin <2
	DT = round(10*(ic.T_s_0-ic.T_a));
	RH = round(ic.RH);
	filename= sprintf('timehistory_%dK_%d_RH',DT,RH);
end
save(sprintf('%s.mat',filename),'timehistory')
cd(my_dir)


























