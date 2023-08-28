#!/bin/sh
#SBATCH --mem=110Gb
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --qos=preemptable
#SBATCH --array=1-45

idx=${SLURM_ARRAY_TASK_ID}
ml R
Rscript /projects/yoki5348/simul_results/perform_ldpred_all.R $idx 


