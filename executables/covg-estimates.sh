#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
#SBATCH --job-name=covg-estimates
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# coverage estimates with jgi_summarize_contigs from metabat2

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
ilm_read_path="${project_path}/00_illumina_for_polishing/illumina_qced"
np_read_path="${project_path}/03_porechop_trimmed/AD_SIP_trimmed_combined.fastq"
assembly_path="${project_path}/re_binning/final_np_polished_contigs_modified.fasta"
out_path="${project_path}/re_binning"

# modules
module load StdEnv/2020 gcc/9.3.0 bowtie2 samtools minimap2 metabat/2.14

# get coverage for each sample with jgi_summarize_contigs
jgi_summarize_bam_contig_depths --outputDepth ${out_path}/r2-np-assembly-samples-depth.txt ${out_path}/mappingResults/*.sorted.bam
