#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=168:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=nuc_div
#SBATCH --error=logs/nuc_div.err.txt
#SBATCH --output=logs/nuc_div.out.txt
#SBATCH --export=All
#SBATCH --array=1

module load angsd/0.935-GCC-10.2.0
sfs_dir=/lustre/home/do271/Trout/BGI/angsd/results/angsd/no_maf_filter/safs
REALSFS=realSFS

for POP in AGA AGB AVB AVD CAM CWA CWB CWC HAY1 NEV2 OWE TYA TYD TYE WEA WEF YHA YHE
do
        echo $POP
        realSFS saf2theta $POP.saf.idx -sfs $POP.sfs -outname $POP
done

THETASTAT=thetaStat

for POP in AGA AGB AVB AVD CAM CWA CWB CWC HAY1 NEV2 OWE TYA TYD TYE WEA WEF YHA YHE
do
	$THETASTAT print $POP.thetas.idx 2>/dev/null | head 
	$THETASTAT do_stat $POP.thetas.idx -win 50000 -step 50000 -outnames $POP.thetas.windows.gz
done
