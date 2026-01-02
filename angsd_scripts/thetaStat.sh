#!/bin/bash

#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=48:00:00
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=saf2theta
#SBATCH --error=logs/theta.err.txt
#SBATCH --output=logs/theta.out.txt
#SBATCH --export=All
#SBATCH --array=1

module load angsd/0.935-GCC-10.2.0
sfs_dir=/lustre/home/do271/Trout/BGI/angsd/results/angsd/no_maf_filter/safs_2

module load angsd/0.935-GCC-10.2.0

### script taken from here and adapted: https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial4_summary_stats/markdowns/summary_stats.md#nucleotide-diversity

### get thetas from sfs 

cd /lustre/home/do271/Trout/BGI/angsd/results/angsd/no_maf_filter/safs_2

for POP in AGA AGB AVB AVD CAM CWA CWB CWC HAY1 NEV2 OWE TYA TYD TYE WEA WEF YHA YHE 
	do realSFS saf2theta $POP.saf.idx -sfs $POP.fold.sfs -outname $POP 
	done

### do windowed thetaStats 
### go for 10k windows - slightly better resolution for plotting whats happening in windows 

for POP in AGA AGB AVB AVD CAM CWA CWB CWC HAY1 NEV2 OWE TYA TYD TYE WEA WEF YHA YHE
 	do thetaStat do_stat $POP.thetas.idx -win 10000 -step 10000 -outnames $POP.thetas.windows.gz
	done
