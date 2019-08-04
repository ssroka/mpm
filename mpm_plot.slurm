#!/bin/bash
#
# Job 
#SBATCH --job-name=mpm_plot
# 

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00
#
#SBATCH --partition=sched_mit_hill
# #SBATCH --mem-per-cpu=50
#SBATCH --constraint=centos7
#SBATCH --time=02:00:00
#SBATCH --tasks-per-node=1
#SBATCH --nodes=1

# Streams
#SBATCH -o mpm_plot_%j.out
#SBATCH -e mpm_plot_%j.err

# Modules
module add mit/matlab/2017b

cd /home/ssroka/sandbox/
# Content
matlab -nodisplay -r "new_calc_CK_CD_driver"












































