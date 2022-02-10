#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=10:0:0
#SBATCH --job-name=11-interproscan
#SBATCH --output=%x.out
#SBATCH --mail-user=ziels@mail.ubc.ca
#SBATCH --mail-type=ALL


#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
assembly_path="${project_path}/07_racon_ilm_polished/"
gene_calls="${project_path}/08_gene-calls"
out_path="${project_path}/11_annotation"
asm="raconx1_medakax1_raconilmx1"


#prepare environment
module load interproscan
mkdir ${out_path}
mkdir ${out_path}/${asm}
mkdir ${out_path}/${asm}/interproscan

## map annotate assembly

interproscan.sh  -cpu 16 -f TSV -f GFF3 -i ${gene_calls}/${asm}/${asm}_proteins.faa -d ${out_path}/${asm}/interproscan

