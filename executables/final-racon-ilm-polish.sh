#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=10:0:0
#SBATCH --job-name=final-racon-ilm-polish
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL


#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
read_path="${project_path}/00_illumina_for_polishing/illumina_qced/racon_ilm_input/R2Sept2020_qced.renamed.interleaved.fastq"
assembly_path="${project_path}/re_polishing/medaka_x3"
out_path="${project_path}/re_polishing/racon_ilm_polishing"

#prepare environment
module load racon bowtie2

## map illumina reads
bowtie2-build ${assembly_path}/consensus.fasta ${assembly_path}/consensus
bowtie2 -x ${assembly_path}/consensus -U ${read_path} -q --no-unal --sensitive -p 16 -S ${out_path}/ilm_to_medaka_consensus_x1.sam


## run racon
racon -t 16 -u \
${read_path} \
${out_path}/ilm_to_medaka_consensus_x1.sam \
${assembly_path}/consensus.fasta > ${out_path}/contigs_medaka.consensus_racon.ilm_x1.fasta