#!/bin/sh
#SBATCH --mem=110gb
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --job-name=simul_gc
#SBATCH --qos=preemptable

ml R
ml singularity
idx=$1

Rscript /pl/active/KellerLab/Yongkang/AM_source_code/perform_ldpred_internal.R $idx


idx=38
for idx in {28,29,21,46,11,15,23,20,33,42,43,44,31,18,41,22,32}
do
sbatch /pl/active/KellerLab/Yongkang/AM_source_code/perform_ldpred_internal.sh $idx
done
