#!/bin/sh
#SBATCH --mem=300Gb
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --qos=preemptable

chr=$1

ml R
for idx in {1..50}
do
Rscript /pl/active/KellerLab/Yongkang/AM_source_code/make_hapmap3_hapfile.R $chr $idx
done
