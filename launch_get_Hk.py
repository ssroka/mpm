#!/bin/python

import subprocess
import os
import fileinput

# the location of this script
base_dir = os.path.dirname(os.path.realpath(__file__))

tof_list = range(20)

for i in range(len(tof_list)):	
	shscript = subprocess.Popen(["sbatch","get_Hk.slurm",str(tof_list[i])],stdin=subprocess.PIPE)
	shscript.wait()








