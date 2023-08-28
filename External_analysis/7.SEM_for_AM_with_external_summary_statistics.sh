#!/bin/sh
#SBATCH --mem=110gb
#SBATCH --array=6-47
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --job-name=simul_gc
#SBATCH --qos=preemptable


ml singularity
idx=${SLURM_ARRAY_TASK_ID}

singularity run /projects/lessem/singularity/openmx.sif Rscript /pl/active/KellerLab/Yongkang/AM_source_code/SEM_for_AM_with_external_summary_statistics.R $idx


