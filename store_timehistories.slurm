#!/bin/bash
#
# Job 
#SBATCH --job-name=store_timehistories 

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00
#
# SBATCH --partition=sched_mit_hill
#SBATCH --partition=sched_mit_twcronin
# #SBATCH --mem-per-cpu=50
#SBATCH --constraint=centos7
#SBATCH --time=12:00:00
#SBATCH --tasks-per-node=1
#SBATCH --nodes=1

# Streams
#SBATCH -o store_timehistories_%j.out
#SBATCH -e store_timehistories_%j.err

# Modules
module add mit/matlab/2017b

cd /home/ssroka/mpm/sandbox/util
# Content
matlab -nodisplay -r "batch_store_timehistories"












































