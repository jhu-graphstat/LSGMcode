#!/bin/bash
#SBATCH -n 1
#SBATCH -o /n/home15/dsussman/log/out-%a.txt
#SBATCH -e /n/home15/dsussman/log/err-%a.txt
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu=3000
#SBATCH -t 600
#SBATCH -a 1-40

module load math/matlab-R2014b
matlab -nojvm -nodisplay -nodesktop -r \
	"job=$SLURM_ARRAY_TASK_ID;addpath(genpath('.'));testTwitter;quit"