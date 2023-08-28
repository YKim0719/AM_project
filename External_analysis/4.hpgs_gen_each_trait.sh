#!/bin/sh
#SBATCH --mem=50Gb
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --qos=preemptable
#SBATCH --array=1-40

idx=${SLURM_ARRAY_TASK_ID}

ml anaconda
conda activate mycustomenv

Rscript /projects/yoki5348/hpgs_gen_each_trait_ver2.R $idx


