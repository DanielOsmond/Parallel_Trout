#!/bin/bash

#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=48:00:00
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=fst
#SBATCH --error=logs/fst.err.txt
#SBATCH --output=logs/fst.out.txt
#SBATCH --export=All
#SBATCH --array=1

module load angsd/0.935-GCC-10.2.0

## This script is used to get pairwise Fst estimates from angsd for each population / group pair

SAFDIR=/lustre/home/do271/Trout/BGI/angsd/results/angsd/no_maf_filter/safs_2 #  Path to per population saf.gz files. An example for the Greenland cod data is: /workdir/cod/greenland-cod/angsd/popminind2/
SAMPLETABLE=/lustre/home/do271/Trout/BGI/angsd/sample_table_3.tsv # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se. An example of such a sample table is: /workdir/cod/greenland-cod/sample_lists/sample_table.tsv
POPCOLUMN=5 # The column index of the variable that you want to group by in the sample table above. In the Greenland project, it's the fifth column, and thus 5
#BASENAME=$4 # Base name of the saf files excluding ".saf.gz". It will be used as the base name of all output files. An example from the Greenland cod project is _global_snp_list_bam_list_realigned_mindp161_maxdp768_minind97_minq20_popminind2
#REALSFS=${5:-/workdir/programs/angsd0.935/angsd/misc/realSFS} # Path to the realSFS module in ANGSD. Default value is /workdir/programs/angsd0.935/angsd/misc/realSFS
THREADS=16 # Number of parallel threads to use, default value is 8.
EXTRA_ARG=${7:-''} # Extra arguments for the SFS estimation step, default value is ''
OUT=/lustre/home/do271/Trout/BGI/angsd/results/angsd/no_maf_filter/fst
cd $SAFDIR

REALSFS=realSFS

I=1
#for POP1 in `tail -n +2 $SAMPLETABLE | cut -f $POPCOLUMN | sort | uniq`; do 
for POP1 in AGA ; do
	J=1
	for POP2 in `tail -n +2 $SAMPLETABLE | cut -f $POPCOLUMN | sort | uniq`; do 
		if [ $I -lt $J ]; then
			echo $POP1'_'$POP2
			if [ ! -f $POP1'.saf.idx' ] || [ ! -f $POP2'.saf.idx' ]; then
				echo 'One or both of the saf.idx files do not exist. Will proceed to the next population pair.'
			else
				# Check if Fst output already exists
				if [ ! -f $POP1'_'$POP2'.fst' ]; then
					# Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)				
					$REALSFS $POP1'.saf.idx' $POP2'.saf.idx' -P $THREADS  > $OUT/$POP1'_'$POP2'.2dSFS'
					# Estimating Fst in angsd
					$REALSFS fst index  $POP1'.saf.idx' $POP2'.saf.idx' -sfs $OUT/$POP1'_'$POP2'.2dSFS' -fstout $OUT/$POP1'_'$POP2'.alpha_beta'
					$REALSFS fst print $OUT/$POP1'_'$POP2'.alpha_beta.fst.idx' > $OUT/$POP1'_'$POP2'.alpha_beta.txt'
					awk '{ print $0 "\t" $3 / $4 }' $OUT/$POP1'_'$POP2'.alpha_beta.txt' > $OUT/$POP1'_'$POP2'.fst'
				fi
				# Check if average Fst output already exists
				if [ ! -f $OUT/$POP1'_'$POP2'.average_fst.txt' ]; then
					# Estimating average Fst in angsd
					$REALSFS fst stats $OUT/$POP1'_'$POP2'.alpha_beta.fst.idx' > $OUT/$POP1'_'$POP2'.average_fst.txt' 
				fi
			fi
		fi
		J=$((J+1))
	done
	I=$((I+1))
done
