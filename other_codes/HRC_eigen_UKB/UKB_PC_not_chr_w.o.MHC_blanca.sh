#!/bin/sh
#SBATCH --mem=180gb
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=simul_gc
#SBATCH --qos=preemptable

chr=$1

~/plink2 --bfile /pl/active/KellerLab/Yongkang/orig_UKB_exclude_HMC_region --threads 20 --pca approx 30 --not-chr "$chr" --out /rc_scratch/yoki5348/UK_PC_each_chr/not_chr"$chr"

for chr in {1..22}
do
sbatch /projects/yoki5348/simul_results/UKB_PC_not_chr_w.o.MHC_blanca.sh $chr
done


###############Remove more variants near in +-2Mb
##############Internal analysis with leave-one chromosome results
###############