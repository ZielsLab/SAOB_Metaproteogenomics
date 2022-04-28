#!/bin/bash

#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
#SBATCH --job-name=bowtie2-build-index
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
fasta="${project_path}/re_binning/combined_bin_set/mappingResults/bins/all_SAOB_bins.fasta"
out_path="${project_path}/re_binning/combined_bin_set/mappingResults/bins/bt2"

# load module
module load bowtie2 

# index command
bowtie2-build ${fasta} ${out_path}/all_SAOB_bins.fasta