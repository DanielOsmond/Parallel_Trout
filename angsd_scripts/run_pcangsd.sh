#!/bin/bash

#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=168:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=angsd_pca
#SBATCH --error=logs/angsd_pca.err.txt
#SBATCH --output=logs/angsd_pca.out.txt
#SBATCH --export=All
#SBATCH --array=1

module purge
module load angsd/0.921-foss-2018a
source ~/.bashrc

master=/lustre/home/do271/Trout/BGI

## This script is used to run PCA using pcangsd. It can be used to run individual-based PCA, estimate selection, inbreeding coefficient, kinship, admixture, and others. The input is a beagle formatted genotype likelihood file.
## See https://github.com/Rosemeis/pcangsd for details

BASEDIR=$master/angsd/results/pcangsd # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
BEAGLE=/lustre/home/do271/Trout/BGI/angsd/results/angsd/pcangsd/merged_filtered.beagle.gz # Path to the beagle formatted genotype likelihood file
MINMAF=0.05 # Minimum allele frequency filter
ANALYSIS=admix # Type of analysis with pcangsd. It can be one of the following: pca, selection, inbreedSites, kinship, admix
MINE=1 # Minimum number of eigenvectors to use in the modelling of individual allele frequencies (relevant for admix)
MAXE=17 # Maximum number of eigenvectors to use in the modelling of individual allele frequencies (relevant for admix)

PREFIX=unlinked_admix

if [ $ANALYSIS = pca ]; then
	pcangsd --beagle $BEAGLE --minMaf $MINMAF --threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX

elif [ $ANALYSIS = selection ]; then 
	pcangsd --beagle $BEAGLE --selection --minMaf $MINMAF --threads 16 -o $BASEDIR'angsd/pcangsd/pcangsd_'$PREFIX --sites_save

elif [ $ANALYSIS = inbreedSites ]; then 
	python2 /workdir/programs/pcangsd/pcangsd.py -beagle $BEAGLE -inbreedSites -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX -sites_save

elif [ $ANALYSIS = kinship ]; then 
	python2 /workdir/programs/pcangsd/pcangsd.py -beagle $BEAGLE -kinship -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX

elif [ $ANALYSIS = admix ]; then 
	for ((E = $MINE; E <= $MAXE; E++)); do 
	echo $E;
	pcangsd -b $BEAGLE --admix -e $E --minMaf $MINMAF --threads 16 -o $BASEDIR'angsd_'$PREFIX'_e'$E;         
	done
fi
