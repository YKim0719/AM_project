#!/bin/sh
#SBATCH --mem=110Gb
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --time=4:00:00

chr=$1; 
pheno_idx=$2
ml R

Rscript /projects/yoki5348/perform_regenie_train_test.R $chr $pheno_idx

for pheno_idx in {1,2,3,4,5,14,18,19,20}
for pheno_idx in {45,38,11,13,12,17,21,22,23,24,25,26,27,29,30}
do
for chr in {1..22}
do
sbatch /projects/yoki5348/perform_regenie_train_test.sh $chr $pheno_idx
    squeue -u yoki5348 > /pl/active/KellerLab/Yongkang/GeneEvolve/queue2.txt  
    filename2='/pl/active/KellerLab/Yongkang/GeneEvolve/queue2.txt'
    nb_lines2=$(cat ${filename2} | wc -l)
      while  (($nb_lines2 >=60)) ;
      do
      sleep 5s
      squeue -u yoki5348 > /pl/active/KellerLab/Yongkang/GeneEvolve/queue2.txt
      filename2='/pl/active/KellerLab/Yongkang/GeneEvolve/queue2.txt'
      nb_lines2=$(cat ${filename2} | wc -l)
      done
done
rm slurm*
rm core*
done
done

for pheno_idx in {31,32,33,34,35,36,37,39,40,41,42,43,44,15,7,8,9}
do
for chr in {1..22}
do
sbatch /projects/yoki5348/perform_regenie_train_test_preemptable.sh $chr $pheno_idx
    squeue -u yoki5348 > /pl/active/KellerLab/Yongkang/GeneEvolve/queue3.txt  
    filename2='/pl/active/KellerLab/Yongkang/GeneEvolve/queue3.txt'
    nb_lines2=$(cat ${filename2} | wc -l)
      while  (($nb_lines2 >=60)) ;
      do
      sleep 5s
      squeue -u yoki5348 > /pl/active/KellerLab/Yongkang/GeneEvolve/queue3.txt
      filename2='/pl/active/KellerLab/Yongkang/GeneEvolve/queue3.txt'
      nb_lines2=$(cat ${filename2} | wc -l)
      done
done
rm slurm*
rm core*
done
done
