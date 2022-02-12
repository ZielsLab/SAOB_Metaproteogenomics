#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
#SBATCH --job-name=09-mapping-exp-sample
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
read_path="${project_path}/00_illumina_for_polishing/illumina_qced/"
assembly_path="${project_path}/07_racon_ilm_polished/"
bowtie_path="${project_path}/09_mapping"
out_path="${project_path}/09_exp_mapping"
asm="raconx1_medakax1_raconilmx1"

# modules
module load bowtie2 samtools

# map illumina reads
bowtie2 -x ${bowtie_path}/${asm}/contigs_index -1 ${read_path}/R2Sept2020_R1.qced.fastq -2 ${read_path}/R2Sept2020_R2.qced.fastq -q --very-sensitive -p 16 -S ${out_path}/${asm}/R2Sept2020.sam

# samtools view, sort, index
samtools view -@ 12 -S -b ${out_path}/${asm}/R2Sept2020.sam > ${out_path}/${asm}/R2Sept2020.bam
samtools sort -@ 12 ${out_path}/${asm}/R2Sept2020.bam -o ${out_path}/${asm}/R2Sept2020}.sorted.bam
samtools index ${out_path}/${asm}/R2Sept2020.sorted.bam
