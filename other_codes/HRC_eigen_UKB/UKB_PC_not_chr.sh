#!/bin/sh
#SBATCH --mem=200gb
#SBATCH --ntasks=60
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --job-name=simul_gc
#SBATCH --qos=normal

chr=$1

~/plink2 --bfile /pl/active/IBG/UKBiobank/GENO/QCed/genotyped/white/eur.qc --threads 20 --pca approx 30 --not-chr "$chr" --out /pl/active/KellerLab/Yongkang/UK_PC_each_chr/not_chr"$chr"_ver2

