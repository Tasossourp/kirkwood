#!/bin/bash

#SBATCH -J kirkwood          # Job name
#SBATCH -o kirkwood.out   # Specify stdout output file (%j expands to jobId)
#SBATCH -p bigmem              # Partition/Queue name
#SBATCH -C skylake               # select 'skylake' architecture
#SBATCH -N 1                     # Total number of nodes requested (32 cores/node)
#SBATCH -t 90:00:00              # Run time (hh:mm:ss) - 0.5 hours
 
#SBATCH -A m2_trr146           # Specify allocation to charge against

#SBATCH --array=0-19
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mail-user=asourpis@uni-mainz.de
#SBATCH --mail-type=ALL

module load bio/GROMACS/2018.1-intel-2018.02

declare -a radius

radius=(0.01, 0.26, 0.51, 0.76, 1.01, 1.26, 1.51, 1.76, 2.01, 2.26, 2.51, 2.76, 3.01, 3.26, 3.51, 3.76, 4.01, 4.26, 4.51, 4.76)

export r=${radius[${SLURM_ARRAY_TASK_ID}]}

srun index.sh &&
srun trj_ch2o_clccn.sh && 
srun trj_ch2o_clh2o.sh
#python3 kirkwood.py

echo 'done'
