#!/bin/bash

#SBATCH -D .
#SBATCH -p sq
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=filt_bam
#SBATCH --error=logs/filt_bam.err.txt
#SBATCH --output=logs/filt_bam.out.txt
#SBATCH --export=All
#SBATCH --array=1

module purge
module load SAMtools/1.9-foss-2018b 

bamlist=/lustre/home/do271/Trout/BGI/freebayes_RG/merged_bams/bam_list.txt
out=/lustre/home/do271/Trout/BGI/freebayes_RG/merged_bams/chr25

# Loop through each BAM file in the list
while IFS= read -r bam_file
do
    # Extract the filename (without extension) from the path
    filename=$(basename "$bam_file")
    filename_noext="${filename%.*}"

    # Use samtools to filter the BAM file and name the output as requested
    samtools view -b "$bam_file" NC_042981.1:23000000-26000000 > $out/"${filename_noext}_chr25:23-26.bam"

done < $bamlist

cd $out
samtools index *
