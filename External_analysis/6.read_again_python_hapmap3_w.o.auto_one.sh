#!/bin/sh
#SBATCH --mem=110Gb
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --qos=preemptable


chr=$1
idx=$2
ml anaconda
conda activate python3conda

python /projects/yoki5348/read_again_by_python_hapmap3_bulk.py $chr $idx


