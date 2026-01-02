#!/bin/bash
#SBATCH -D .
#SBATCH -p highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=angsd_snps
#SBATCH --error=logs/GL.err.txt
#SBATCH --output=logs/GL.out.txt
#SBATCH --export=All
#SBATCH --array=1

module purge
module load angsd/0.921-foss-2018a
module load SAMtools/1.7-foss-2018a

master=/lustre/home/do271/Trout/BGI
DATA=$master/freebayes_RG/merged_bams
REF=/gpfs/ts0/home/do271/Trout/reference/GCF_901001165.1_fSalTru1.1_genomic.fna

### This command calculates the Genotype likelihoods ###
### Use the whole dataset ###
### Applies mapping quality ###
### Output in binary -doGlf 1 (change this option to -doGlf 4 for txt output ###
### -nThreads - me trying to speed it up, dunno if this works ###  

angsd -b  $DATA/bam_list.txt -ref $REF -out $master/angsd/results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 5 -doCounts 1 \
	-GL 2 -doGlf 1\
	-nThreads 32

### Command to call genotypes ### 
### doMajorMinor 1 - Infer major and minor from GL ###
### doPost 2 - Uniform prior. Can set as 1 to use frequency as prior (relies on HWE, not true for sep pops) ###


angsd -glf $master/angsd/results/ALL.glf.gz -fai $REF.fai -nInd 137 -out $master/angsd/results/ALL \
	-doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1
