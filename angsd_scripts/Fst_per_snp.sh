#!/bin/bash

#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=snp_fst
#SBATCH --error=logs/snp_fst.err.txt
#SBATCH --output=logs/snp_fst.out.txt
#SBATCH --export=All
#SBATCH --array=1

module load angsd/0.935-GCC-10.2.0

pop_pairs=/lustre/home/do271/Trout/BGI/angsd/results/sample_lists/pop_pairs.txt
saf_dir=/lustre/home/do271/Trout/BGI/angsd/results/angsd/safs/results
fst_dir=$saf_dir/per_locus
REALSFS=realSFS

#simpleID_array=( `cat $pop_pairs | cut -f 1` )
#POP1=${simpleID_array[(($SLURM_ARRAY_TASK_ID))]}

#pop2_array=( `cat $pop_pairs | cut -f 2` )
#POP2=${pop2_array[(($SLURM_ARRAY_TASK_ID))]}

POP1=AGA
POP2=AGB

$REALSFS $saf_dir/$POP1.saf.idx $saf_dir/$POP2.saf.idx -P 16 > $saf_dir/${POP1}_${POP2}.2dSFS

$REALSFS fst index  $saf_dir/$POP1.saf.idx $saf_dir/$POP2.saf.idx -sfs $saf_dir/${POP1}_${POP2}.2dSFS -fstout \
$fst_dir/${POP1}_${POP2}.alpha_beta

$REALSFS fst print $fst_dir/${POP1}_${POP2}.alpha_beta.fst.idx > ${POP1}_${POP2}.alpha_beta.txt
