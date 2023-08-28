#!/bin/sh
#SBATCH --mem=200gb
#SBATCH --ntasks=60
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --job-name=simul_gc
#SBATCH --qos=normal

chr=$1

~/plink2 --bfile /pl/active/KellerLab/Yongkang/orig_UKB_exclude_HMC_region --threads 20 --pca approx 30 --not-chr "$chr" --out /scratch/alpine/yoki5348/UK_PC_each_chr/not_chr"$chr"

for chr in {1..22}
do
sbatch /projects/yoki5348/simul_results/UKB_PC_not_chr_w.o.MHC.sh $chr
done

