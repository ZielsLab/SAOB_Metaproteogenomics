#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --array=1-15
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:0:0
#SBATCH --job-name=finalBins_mapping.sh
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# array jobs 
mapping_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" saob_metagenomes.txt)

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
ref_path="${project_path}/re_binning/combined_bin_set/mappingResults/bins/bt2/all_SAOB_bins.fasta"
reads_path=$(dirname $mapping_file)
sample_name=$(basename $mapping_file _R1.qced.fastq)
r1_file="${reads_path}/${sample_name}_R1.qced.fastq"
r2_file="${reads_path}/${sample_name}_R2.qced.fastq"
out_path="${project_path}/re_binning/combined_bin_set/mappingResults"
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
