#!/bin/sh
#SBATCH --mem=40gb
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --job-name=simul_gc
#SBATCH --qos=normal

idx=$1
ml anaconda
conda activate mycustomenv
Rscript /projects/yoki5348/simul_results/UKB_PC_adj_each_chr.R $idx



for idx in {7..47}
do
sbatch /projects/yoki5348/simul_results/UKB_PC_adj_each_chr.sh $idx
done

for chr in {1..22}
do
sbatch /projects/yoki5348/simul_results/simul_PC_gen1_not_chr.sh $chr
done
