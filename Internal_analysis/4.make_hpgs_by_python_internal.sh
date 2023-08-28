#!/bin/sh
#SBATCH --mem=110Gb
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --qos=preemptable
#SBATCH --array=1-50

chr=$1
dirname=$2
idx=${SLURM_ARRAY_TASK_ID}
ml anaconda
conda activate python3conda


python /projects/yoki5348/make_hpgs_by_python_internal_ver2.py $chr $idx $dirname

