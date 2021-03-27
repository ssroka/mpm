#!/bin/python

import subprocess
import os
import fileinput

# the location of this script
base_dir = os.path.dirname(os.path.realpath(__file__))

# RH_list = [88.0, 90.0, 92.0, 94.0, 96.0, 98.0]
RH_list = [80.0, 82.0, 84.0, 86.0]
RH_list = [80.0, 82.0, 84.0, 86.0]
DT_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
Ts_list = [27.0, 28.0, 29.0]

#for i_RH in range(len(RH_list)):	
for i_RH in [2, 3]:	
	for i_DT in range(len(DT_list)):	
		for i_Ts in range(len(Ts_list)):	
			#exp_num = 100*(2+i_Ts) + 10*(i_DT) + 1*(1+i_RH)
			exp_num = 1000*(5+i_Ts) + 100*(i_DT) + 10*(1+i_RH)
			shscript = subprocess.Popen(["sbatch","mpm.slurm",str(exp_num)],stdin=subprocess.PIPE)
			shscript.wait()








