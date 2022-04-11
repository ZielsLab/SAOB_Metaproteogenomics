#!/bin/bash
#SBATCH --account=def-shallam
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=6:0:0
#SBATCH --job-name=kofam-annotate.sh
#SBATCH --output=kofam-annotate.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# paths 
project_path="/home/eamcdani/projects/def-shallam/eamcdani/sakinaw"
bins_path="${project_path}/final_bins/pathway_annotation"
kofam_path="/home/eamcdani/bin/kofam_scan-1.3.0"

# load modules 
module load ruby/2.7.1 hmmer 

# cd to directory
cd ${bins_path}

# kofam koala for loop for each genome FAA file
for fasta in *.faa; do
	name=$(basename $fasta .faa);
	${kofam_path}/exec_annotation $fasta -o $name-kofam-annotations.txt -p ${kofam_path}/profiles/ -k ${kofam_path}/ko_list --cpu 8;
done 