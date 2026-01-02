#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=168:00:00
#SBATCH --nodes=25
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=beag
#SBATCH --error=logs/get_beagle_glf4.err.txt
#SBATCH --output=logs/get_beagle_glf4.out.txt
#SBATCH --export=All
#SBATCH --array=1

module purge
module load angsd/0.921-foss-2018a

master=/lustre/home/do271/Trout/BGI

BAMLIST=$master/freebayes_RG/merged_bams/bam_list.txt # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
BASEDIR=$master/angsd/results/ # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
REFERENCE=/gpfs/ts0/home/do271/Trout/reference/GCF_901001165.1_fSalTru1.1_genomic.fna # Path to reference genome
SNPLIST=$master/angsd/results/angsd/snp_lists/high_conf_snps.txt  #path to the SNP list

## Build base name of output files
OUTBASE=high_conf_snps
CHRLIST=high_conf_snps.chrs

#angsd sites index /lustre/home/do271/Trout/BGI/angsd/results/angsd/not_problem_snps_2.txt

## Get beagle formatted genotype likelihood
angsd \
-b $BAMLIST \
-anc $REFERENCE \
-out /lustre/archive/scratch/do271/angsd/$OUTBASE \
-GL 1 \
-dovcf 1 \
-doGlf 4 \
-doMaf 1 \
-doPost 2 \
-doMajorMinor 3 \
-P 32 \
-sites $SNPLIST \
-rf $master/angsd/results/angsd/snp_lists/high_conf_snps.chrs \
>& $BASEDIR'nohups/'$OUTBASE'_get_beagle_2.log'
