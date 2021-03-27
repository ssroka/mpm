#!/bin/bash
#
#SBATCH -J mm
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o mm.out
#SBATCH -e mm.err
# SBATCH --partition=sched_mit_hill
#SBATCH --partition=sched_mit_twcronin
 
. /etc/profile.d/modules.sh
module add mit/matlab/2016a
cd /home/ssroka/sandbox
matlab -nodesktop -nosplash -r "dropletEvolution_driver($SLURM_JOB_ID)"












