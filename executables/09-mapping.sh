#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
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
mkdir ${out_path}
mkdir ${out_path}/${asm}

#anvio rename contigs
module load python scipy-stack diamond hmmer prodigal bowtie2
source ~/projects/rrg-ziels/shared_tools/virtual_envs/anvio_7.0/bin/activate

anvi-script-reformat-fasta ${assembly_path}/contigs_medaka.consensus_racon.ilm_x1.fasta -o ${assembly_path}/contigs_medaka.consensus_racon.ilm_x1.renamed.fasta -l 0 --simplify-names

deactivate

## map illumina reads
bowtie2-build ${assembly_path}/contigs_medaka.consensus_racon.ilm_x1.renamed.fasta ${out_path}/${asm}/contigs_index

for dir in ${read_path}/*_R1.qced.fastq; do
        dir=${dir%*_R1.qced.fastq}
        sample=`basename "$dir"`
        echo $sample

bowtie2 -x ${out_path}/${asm}/contigs_index -1 ${read_path}/${sample}_R1.qced.fastq -2 ${read_path}/${sample}_R2.qced.fastq -q --very-sensitive -p 16 -S ${out_path}/${asm}/${sample}.sam

samtools view -b -S ${out_path}/${asm}/${sample}.sam > ${out_path}/${asm}/${sample}.bam


done

