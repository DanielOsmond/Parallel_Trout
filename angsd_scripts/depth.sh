#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=12:00:00
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T111858
#SBATCH --job-name=depth
#SBATCH --error=logs/depth.err.txt
#SBATCH --output=logs/depth.out.txt
#SBATCH --export=All
#SBATCH --array=1

module purge
module load SAMtools/1.9-foss-2018b 

sample_table=/lustre/home/do271/Trout/BGI/angsd/sample_table_3.tsv
unique_prefixes=$(awk '{ if (NR > 1) print $1 }' $sample_table | sort | uniq)
#bamlists=/lustre/home/do271/Trout/BGI/angsd/results/sample_lists/bam_list_per_pop

cd /lustre/home/do271/Trout/BGI/freebayes_RG/merged_bams/chr25

# Loop through each unique prefix
#for prefix in $unique_prefixes
#do
#samtools depth -a -b golden.bed ${prefix}.merged.dups_chr25:23-26.bam > ${prefix}.chr25.depth.txt
#done

#for pop in AGA AGB AVB AVD CAM CWA CWB CWC HAY1 NEV OWE TYA TYD TYE WEA WEF YHA YHE
#do samtools depth -b golden.bed -f ${pop}_bam_list.txt  > depth/$pop.chr25.4x.depth.txt
#done

for pop in HAY1 NEV2
do samtools depth -b golden.bed -f ${pop}_bam_list.txt  > depth/$pop.chr25.4x.depth.txt
done

