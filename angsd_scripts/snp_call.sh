#!/bin/bash
#SBATCH -D .
#SBATCH -p highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=angsd_snps
#SBATCH --error=logs/snp_call.err.txt
#SBATCH --output=logs/snp_call.out.txt
#SBATCH --export=All
#SBATCH --array=1

module purge
module load angsd/0.921-foss-2018a
module load SAMtools/1.7-foss-2018a

master=/lustre/home/do271/Trout/BGI
DATA=$master/freebayes_RG/merged_bams
REF=/gpfs/ts0/home/do271/Trout/reference/GCF_901001165.1_fSalTru1.1_genomic.fna

angsd -b  $DATA/bam_list.txt -ref $REF -out $master/angsd/results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 5 -doCounts 1 \
	-GL 2 -doGlf \
	-nThreads 32
