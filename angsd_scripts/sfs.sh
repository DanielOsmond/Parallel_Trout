#!/bin/bash

#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=realSFS
#SBATCH --error=logs/realSFS.err.txt
#SBATCH --output=logs/realSFS.out.txt
#SBATCH --export=All
#SBATCH --array=1

module load angsd/0.935-GCC-10.2.0
saf_dir=/lustre/home/do271/Trout/BGI/angsd/results/angsd/maf_05
sfs_dir=/lustre/home/do271/Trout/BGI/angsd/results/angsd/safs/results
REALSFS=realSFS

for POP in HAY1 CWC
do
        echo $POP
        $REALSFS $saf_dir/${POP}_maf.saf.idx > $sfs_dir/$POP.sfs
done
