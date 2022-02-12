#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
#SBATCH --job-name=10-metabat-binning
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
mapping_path="${project_path}/09_mapping/raconx1_medakax1_raconilmx1"
out_path="${project_path}/10_binning"
assembly_path="${project_path}/07_racon_ilm_polished"

module load StdEnv/2020 gcc/9.3.0 metabat/2.14

# get depth
jgi_summarize_bam_contig_depths --outputDepth ${out_path}/r2-nanopore-assembly-v1-depth.txt ${mapping_path}/*.sorted.bam

# run metabat 
metabat2 -i ${assembly_path}/contigs_medaka.consensus_racon.ilm_x1.renamed.fasta -a ${out_path}/r2-nanopore-assembly-v1-depth.txt -o ${out_path}/nanopore_v1-bin

