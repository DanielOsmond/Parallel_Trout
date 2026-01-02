#!/bin/bash

#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=24:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=LD_plot
#SBATCH --error=LD_plot.err.txt
#SBATCH --output=LD_plot.out.txt
#SBATCH --export=All
#SBATCH --array=1

module load R/4.0.0-foss-2020a

CHR=NC_042981.1
START_POS=20000000
END_POS=30000000

cd /lustre/home/do271/Trout/BGI/angsd/results/ld

cat ld/81.ld | bash LD_blocks.sh 
