#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --array=1-16
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=5:0:0
#SBATCH --job-name=bowtie2-mapping-to-bins
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# array jobs 
mapping_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" saob_metagenomes.txt)

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
ref_path="${project_path}/re_binning/combined_bin_set/all_saob_bins.fasta"
reads_path=
out_path="${project_path}/re_binning/combined_bin_set/mappingResults"
sample_name=
r1_file=
r2_file=
out_name="${out_path}/${sample_name}-vs-bins"

# load modules
module load bowtie2 samtools

# bowtie2 mapping command 
bowtie2 -x $ref_path -1 $r1_file -2 $r2_file -q --very-sensitive -p 8 -S ${out_name}.sam

# samtools SAM to BAM 
samtools view -@ 12 -S -b ${out_name}.sam > ${out_name}.bam

# samtools sort and index
samtools sort -@ 12 ${out_name}.bam -o ${out_name}.sorted.bam
samtools index ${out_name}.sorted.bam
