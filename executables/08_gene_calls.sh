#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=5:0:0
#SBATCH --job-name=08-prodigal.out
#SBATCH --output=%x.out
#SBATCH --mail-user=ziels@mail.ubc.ca
#SBATCH --mail-type=ALL


#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
assembly_path="${project_path}/07_racon_ilm_polished/contigs_medaka.consensus_racon.ilm_x1.fasta"
out_path="${project_path}/08_gene-calls"
asm="raconx1_medakax1_raconilmx1"

#prepare environment
module load prodigal
mkdir ${out_path}
mkdir ${out_path}/${asm}

## run prodigal

prodigal -i ${assembly_path} -o ${out_path}/${asm}/${asm}_coords.gbk -a ${out_path}/${asm}/${asm}_proteins.faa -d ${out_path}/${asm}/${asm}_genes.fasta -p meta
