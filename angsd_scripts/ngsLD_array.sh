#!/bin/bash

#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=168:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=ngsLD
#SBATCH --error=logs/ngsLD_test.err.txt
#SBATCH --output=logs/ngsLD_test.out.txt
#SBATCH --export=All
#SBATCH --array=1-10

module purge 
ngsLD=/lustre/home/do271/software/ngsLD/
module load GCCcore/9.3.0
module load zlib/1.2.11-GCCcore-9.3.0
module load GSL/2.6-GCC-9.3.0

master=/lustre/home/do271/Trout/BGI
BASEDIR=$master/angsd/results

pos=/lustre/home/do271/Trout/BGI/angsd/results/angsd/chrs_snps/pos_list.txt
beagles=/lustre/home/do271/Trout/BGI/angsd/results/angsd/beagle_by_chrs/beagle_list.txt

pos_dir=/lustre/home/do271/Trout/BGI/angsd/results/angsd/chrs_snps
beagle_dir=/lustre/home/do271/Trout/BGI/angsd/results/angsd/beagle_by_chrs

beagleID_array=( `cat $beagles | cut -f 1` )
beagleID=${beagleID_array[(($SLURM_ARRAY_TASK_ID))]} 

pos_ID_array=( `cat $pos | cut -f 1` )
posID=${pos_ID_array[(($SLURM_ARRAY_TASK_ID))]}

n_array=( `cat $pos | cut -f 2` )
n=${n_array[(($SLURM_ARRAY_TASK_ID))]} 

$ngsLD/ngsLD \
--geno $beagle_dir/$beagleID \
--pos $pos_dir/$posID \
--probs \
--n_ind 149 \
--n_sites $n \
--max_kb_dist 0 \
--n_threads 1 \
--out $BASEDIR/angsd/$posID.ld 

