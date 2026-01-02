#!/bin/bash
#SBATCH -D .
#SBATCH -p highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=angsd_snps
#SBATCH --error=logs/angsd_SNP.err.txt
#SBATCH --output=logs/angsd_SNP.out.txt
#SBATCH --export=All
#SBATCH --array=1

module purge
module load angsd/0.921-foss-2018a

master=/lustre/home/do271/Trout/BGI

BAMLIST=$master/freebayes_RG/merged_bams/bam_list.txt # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
BASEDIR=$master/angsd/results/ # Path to the base directory where output files will be written to a subdirectory named "angsd/". An example for the Greenland cod data is: /workdir/cod/greenland-cod/
REFERENCE=/gpfs/ts0/home/do271/Trout/reference/GCF_901001165.1_fSalTru1.1_genomic.fna # Path to reference genome
MINDP=0.3 # Minimum depth filter
MAXDP=15 # Maximum depth filter
MININD=2 # Minimum individual filter
MINQ=20 # Minimum quality filter
MINMAF=0.05 # Minimum minor allele frequency filter
MINMAPQ=20 # Minimum mapping quality (alignment score) filter, default value is 20
ANGSD=angsd # Path to ANGSD, default value is /workdir/programs/angsd0.931/angsd/angsd
THREADS=32 # Number of parallel threads to use, default value is 8.
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50' # Extra arguments when running ANGSD, default value is '-remove_bads 1 -only_proper_pairs 1 -C 50'

## Extract the name of the bam list (excluding path and suffix)
BAMLISTNAME=`echo $BAMLIST | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`

## Build base name of output files
OUTBASE=$BAMLISTNAME'_mindp'$MINDP'_maxdp'$MAXDP'_minind'$MININD'_minq'$MINQ

## Call SNPs
#$ANGSD -b $BAMLIST -ref $REFERENCE -out $BASEDIR'angsd/'$OUTBASE \
#-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 4 -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
#-setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD \
#-minQ $MINQ -minMapQ $MINMAPQ \
#-SNP_pval 1e-6 -minMaf $MINMAF \
#-P $THREADS \
#$EXTRA_ARG \
#>& $BASEDIR'nohups/'$OUTBASE'.log'

## Create a SNP list to use in downstream analyses
gunzip -c $BASEDIR'angsd/'$OUTBASE'.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'
$ANGSD sites index $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'

## Also make it in regions format for downstream analyses
cut -f 1,2 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sed 's/\t/:/g' > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.regions'

## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
cut -f1 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sort | uniq > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.chrs'
