#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=12:0:0
#SBATCH --job-name=09-mapping
#SBATCH --output=%x.out
#SBATCH --mail-user=ziels@mail.ubc.ca
#SBATCH --mail-type=ALL


#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
read_path="${project_path}/illumina_qced"
assembly_path="${project_path}/07_racon_ilm_polished/"
out_path="${project_path}/09_mapping"
asm="raconx1_medakax1_raconilmx1"


#prepare environment
module load samtools

#sam to bam
for dir in ${read_path}/*_R1.qced.fastq; do
        dir=${dir%*_R1.qced.fastq}
        sample=`basename "$dir"`
        echo $sample

samtools sort -@ 4 ${out_path}/${asm}/${sample}.sam -o ${out_path}/${asm}/${sample}.bam
samtools index ${out_path}/${asm}/${sample}.bam 

done

